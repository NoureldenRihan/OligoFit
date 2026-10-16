"""
INSTALLATION:
    pip install primer3-py biopython pandas

USAGE:
    from oligofit import design_oligos, evaluate_oligos_csv

    # Minimal: uses "custom" mode with default parameters
    df, frag_df, primers = design_oligos(sequence="ATGC...")

    # Preset mode: short oligos (40-60bp oligos, 15-25bp overlaps, 400-500bp fragments)
    # Only the three structural ranges are set by the preset; all other params stay user-controllable.
    df, frag_df, primers = design_oligos(
        sequence = "ATGC...",
        mode     = "short_oligos",
    )

    # Preset mode: long oligos (680-700bp oligos, 30-40bp overlaps, 3000-4000bp fragments)
    df, frag_df, primers = design_oligos(
        sequence = "ATGC...",
        mode     = "long_oligos",
    )

    # Full custom parameters (mode="custom" is the default)
    df, frag_df, primers = design_oligos(
        sequence                    = "ATGC...",
        mode                        = "custom",

        # STRUCTURAL
        fragment_length_range       = (400, 500),   # None = single fragment
        oligo_length_range          = (40, 60),     # supports >60bp via sliding window
        overlap_length_range        = (15, 25),
        
        # QUALITY THRESHOLDS
        overlap_tm_range            = (60.0, 65.0),
        tm_warning_buffer           = 1.0,     
        gc_warning_buffer           = 2.0,     
        overlap_gc_range            = (30.0, 65.0),
        max_homopolymer_run         = 4,       # example: AAAA or GGGG
        max_dinucleotide_repeat     = 3,       # example: ATATAT or GCGCGC
        repeat_vicinity_buffer      = 5,       # flag repeats if within ±N bases of overlap
        hairpin_dg_threshold        = -9.0,    # kcal/mol 
        hairpin_dg_buffer           = 1.0,     # kcal/mol 
        homodimer_dg_threshold      = -9.0,    # kcal/mol
        homodimer_dg_buffer         = 1.0,     # kcal/mol
        gc_clamp_length             = 3,
        gc_clamp_min_gc             = 1,
        heterodimer_dg_threshold    = -9.0,    # kcal/mol 
        heterodimer_dg_buffer       = 1.0,     # kcal/mol 

        # MISPRIME DETECTION
        misprime_length             = 20,      # post-tiling probe length
        misprime_tip_length         = 8,       # exact-match probe tip 
        misprime_max_mismatches     = 2,       # allowed mismatches in probe body
        misprime_kmer_size          = 12,      # pre-tiling risk map seed size 
        
        # HETERODIMER DETECTION
        hetdimer_kmer_size          = 12,      # pre-tiling RC risk map seed size
        heterodimer_neighbor_range  = 2,       # ±N neighbors checked (±1 always skipped)
        
        # DP PENALTY BASE VALUES
        penalty_tm                  = 20,
        penalty_gc                  = 8,
        penalty_homopolymer         = 10,
        penalty_dinucleotide        = 6,
        penalty_hairpin             = 10,
        penalty_homodimer           = 6,
        penalty_gc_clamp            = 5,
        penalty_misprime            = 6,
        penalty_hetdimer_map        = 8,       
        penalty_body_homopolymer    = 3,       
        penalty_body_dinucleotide   = 3,       
        penalty_severity_multiplier = 1.0,     

        # FRAGMENT CUT SCORING
        cut_search_window           = 30,
        cut_score_window            = 20,
        cut_score_repeat_penalty    = 1000.0,

        # JUNCTION AND OUTER PRIMERS
        junction_overlap_length     = 25,
        junction_tm_range           = (55.0, 68.0),
        outer_primer_length_range   = (15, 35),
        outer_primer_tm_range       = (55.0, 68.0),

        # ASSEMBLY METHOD
        assembly_method             = "gibson",          # 'stitching_pcr' | 'gibson' | 'golden_gate'
        #   stitching_pcr → fragments share junction_overlap_length bp for overlap-extension PCR
        #   gibson        → junction_overlap ≥40 bp enforced and Tm min 48°C
        #   golden_gate   → Type IIS cassette outer primers and 4-bp overhang design
        # Golden Gate specific options
        goldengate_enzyme           = "BsaI",             # 'BsaI' | 'BsmBI' | 'BbsI' (enzymes auto-switch if internal sites found)
       
        # THERMODYNAMIC CONDITIONS
        mv_conc                     = 50.0,    # mM  Na+
        dv_conc                     = 3.0,     # mM  Mg2+
        dntp_conc                   = 0.8,     # mM
        dna_conc                    = 200.0,   # nM

        # WORKFLOW CONTROL
        run_oligos                  = True,    # False = fragment table only
        run_primers                 = True,    # False = no outer primers
        fragment_indices            = None,    # example: [1, 3] to tile only frags 1 & 3
        n_workers                   = None,    # None = all CPU cores
        output_csv                  = "mySequence.csv",
    )

PRESET MODES:
    "custom"       : (default) all parameters are user-controlled
    "short_oligos" : oligo 40-60bp, overlap 15-25bp, fragment 400-500bp
    "long_oligos"  : oligo 680-700bp, overlap 30-40bp, fragment 3000-4000bp
    Presets lock only structural ranges; all other params remain user-controllable.

OUTPUTS:
    Returns a tuple of three DataFrames:
      df     : one row per PCA oligo (ready to order sequences)
      frag_df: one row per fragment, including junction info to next fragment
      primers: one forward + one reverse outer PCR primer per fragment (Extra Golden Gate primers if chosen as assembly method)

    when using `run_oligos=False`, the first return (df) is an empty DataFrame and frag_df still contains the full fragment table

    Each DataFrame is also saved to a CSV when output_csv is set:
      <output_csv>                  : oligo table
      <output_csv>_fragments.csv    : fragment table
      <output_csv>_outer_primers.csv: outer PCR primer table
"""

import sys
import re
import math
import time
import heapq
import warnings as _warnings
import pandas as pd
import primer3
from Bio.Seq import Seq
from typing import Tuple, Optional, List, Dict, Set
from concurrent.futures import ProcessPoolExecutor, as_completed
import concurrent.futures
from collections import defaultdict
from functools import lru_cache

_RC_TABLE: dict = str.maketrans('ACGTacgt', 'TGCAtgca')

# DEFAULT VALUES

DEFAULT_FRAGMENT_LENGTH_RANGE    = (400, 500)
DEFAULT_OLIGO_LENGTH_RANGE       = (40, 60)
DEFAULT_OVERLAP_LENGTH_RANGE     = (15, 25)
DEFAULT_OVERLAP_TM_RANGE         = (60.0, 65.0)
DEFAULT_TM_WARNING_BUFFER        = 1.0
DEFAULT_GC_WARNING_BUFFER        = 2.0
DEFAULT_HAIRPIN_DG_BUFFER        = 1.0
DEFAULT_HOMODIMER_DG_BUFFER      = 1.0
DEFAULT_HETERODIMER_DG_BUFFER    = 1.0
DEFAULT_OVERLAP_GC_RANGE         = (30.0, 65.0)
DEFAULT_MAX_HOMOPOLYMER_RUN      = 4
DEFAULT_MAX_DINUCLEOTIDE_REPEAT  = 3
DEFAULT_REPEAT_VICINITY_BUFFER   = 5
# Ref: https://eu.idtdna.com/pages/support/faqs/how-do-i-use-the-oligoanalyzer-tool-to-analyze-possible-hairpins-and-dimers-formed-by-my-oligo
DEFAULT_HAIRPIN_DG_THRESHOLD     = -9.0           # kcal/mol

# Thermodynamic reaction conditions match IDT OligoAnalyzer tool defaults (https://www.idtdna.com/calc/analyzer)
# Unit note: IDT displays oligo concentration in µM, primer3 uses nM, so, 0.2 µM = 200 nM
DEFAULT_MV_CONC                  = 50.0           # mM  monovalent salt (Na⁺/K⁺)
DEFAULT_DV_CONC                  = 3.0            # mM  divalent salt (Mg²⁺)
DEFAULT_DNTP_CONC                = 0.8            # mM  dNTPs
DEFAULT_DNA_CONC                 = 200.0          # nM  oligo strand concentration (IDT: 0.2 µM)
DEFAULT_OUTPUT_CSV               = "mySequence.csv"

# Sliding window parameters for structure checks on oligos > 60 bp
MAX_PRIMER3_SEQ_LENGTH = 60   # primer3 Oligo Checks hard cap
SLIDING_WINDOW_SIZE    = 60   # window for long-oligo structure checks
SLIDING_WINDOW_STEP    = 30   # 50% overlap to catch cross-boundary structures
# Oligos longer than this skip internal windows and only check the first + last 60 bp
# (internal hairpins in PCA oligos are transient and resolved during extension)
LONG_OLIGO_TERMINAL_ONLY_THRESHOLD = 200

DEFAULT_MISPRIME_LENGTH          = 20 # window for misprime check
DEFAULT_MISPRIME_TIP_LENGTH      = 8  # exact-match probe tip
DEFAULT_MISPRIME_MAX_MISMATCHES  = 2  # allowed mismatches in probe body
DEFAULT_JUNCTION_OVERLAP_LENGTH  = 25
DEFAULT_JUNCTION_TM_RANGE        = (55.0, 68.0)
ASSEMBLY_METHODS            = ("stitching_pcr", "gibson", "golden_gate")
DEFAULT_ASSEMBLY_METHOD     = "gibson"
GIBSON_MIN_JUNCTION_OVERLAP = 40  
GIBSON_MIN_JUNCTION_TM      = 48.0
GOLDEN_GATE_ENZYMES = {
    "BsaI":  {"recognition": "GGTCTC", "rc_recognition": "GAGACC", "fwd_spacer": "A"},
    "BsmBI": {"recognition": "CGTCTC", "rc_recognition": "GAGACG", "fwd_spacer": "A"},
    "BbsI":  {"recognition": "GAAGAC", "rc_recognition": "GTCTTC", "fwd_spacer": "AA"},
}
DEFAULT_GOLDENGATE_ENZYME         = "BsaI"
GOLDEN_GATE_EXTRA_BASES           = "GCATAC"
GOLDEN_GATE_MIN_OVERHANG_DISTANCE = 2
# Upper bound only. With the minimum overhang distance enforced while choosing cut sites, random native
# sequences were found to support about 13-23 fragments (median ~18) before no valid overhang remains.
MAX_GOLDENGATE_FRAGMENTS          = 50

DEFAULT_HOMODIMER_DG_THRESHOLD   = -9.0   # kcal/mol
DEFAULT_GC_CLAMP_LENGTH          = 3     
DEFAULT_GC_CLAMP_MIN_GC          = 1      
DEFAULT_HETERODIMER_DG_THRESHOLD = -9.0   # kcal/mol
DEFAULT_HETERODIMER_NEIGHBOR_RANGE = 2    # check ±N neighbors (skipping ±1)
DEFAULT_PENALTY_SEVERITY_MULTIPLIER = 1.0
# How Penalty Severity is calculated for each check:
#   Tm out of range    → degrees C outside the range 
#   GC% out of range   → percentage points outside the range 
#   Homopolymer        → extra repetitions beyond max_homopolymer_run
#   Dinucleotide       → extra repetitions beyond max_dinucleotide_repeat
#   Hairpin           → |dG - threshold| in kcal/mol beyond threshold
#   Homodimer       → |dG - threshold| in kcal/mol beyond threshold
#   Heterodimer (pre) → number of off-target hits
#   Misprime (pre) → number of off-target hits

DEFAULT_OUTER_PRIMER_LENGTH_RANGE = (15, 35)
DEFAULT_OUTER_PRIMER_TM_RANGE     = (55.0, 68.0)
DEFAULT_CUT_SEARCH_WINDOW        = 30
MIN_RELIABLE_TM_LENGTH           = 10

# DP penalty weights
DP_PENALTY_TM           = 20   
DP_PENALTY_GC           = 8    
DP_PENALTY_HOMOPOLYMER  = 10   
DP_PENALTY_DINUCLEOTIDE = 6    
DP_PENALTY_HAIRPIN      = 10   
DP_PENALTY_HOMODIMER    = 6    
DP_PENALTY_GC_CLAMP     = 5    
DP_PENALTY_MISPRIME     = 6    
DP_PENALTY_BODY_HOMOPOLYMER  = 3   
DP_PENALTY_BODY_DINUCLEOTIDE = 3
DP_PENALTY_HETDIMER_MAP    = 8    
DP_PENALTY_FORCED       = 1000 
DEFAULT_MISPRIME_RISK_PENALTY = DP_PENALTY_MISPRIME
DEFAULT_MISPRIME_KMER_SIZE    = 12  
DEFAULT_HETDIMER_KMER_SIZE = 12
MISPRIME_SAFE_BUFFER = 5   
HETDIMER_SAFE_BUFFER = 5   

# ── Preset Modes ──────────────────────────────────────────────────
VALID_MODES = ("custom", "short_oligos", "long_oligos")

PRESET_SHORT_OLIGOS = {
    "fragment_length_range": (400, 500),
    "oligo_length_range":    (40, 60),
    "overlap_length_range":  (15, 25),
}

PRESET_LONG_OLIGOS = {
    "fragment_length_range": (3000, 4000),
    "oligo_length_range":    (680, 700),
    "overlap_length_range":  (30, 40),
}

_PRESETS = {
    "short_oligos": PRESET_SHORT_OLIGOS,
    "long_oligos":  PRESET_LONG_OLIGOS,
}


# SECTION 1: BASIC SEQUENCE UTILITIES

@lru_cache(maxsize=None)
def calc_tm(seq: str, mv_conc: float, dv_conc: float,
            dntp_conc: float, dna_conc: float) -> float:
    """
    Calculate the melting temperature (Tm) using SantaLucia 1998 nearest-neighbor via primer3

    Args:
        seq       : DNA sequence (5'→3')
        mv_conc   : monovalent salt in mM (Na+/K+)
        dv_conc   : divalent salt in mM (Mg2+)
        dntp_conc : dNTP concentration in mM
        dna_conc  : oligo strand concentration in nM

    Returns 0.0 for empty sequences and  a UserWarning for sequences < MIN_RELIABLE_TM_LENGTH (model unreliable below this length)
    """
    if len(seq) == 0:
        _warnings.warn("calc_tm() called with empty sequence. Returning 0.0°C.", UserWarning, stacklevel=2)
        return 0.0
    if len(seq) < MIN_RELIABLE_TM_LENGTH:
        _warnings.warn(f"calc_tm(): sequence '{seq}' ({len(seq)} bp) is below the {MIN_RELIABLE_TM_LENGTH} bp reliability threshold.", UserWarning, stacklevel=2)
    return primer3.calc_tm(
        seq.upper(),
        mv_conc=mv_conc,
        dv_conc=dv_conc,
        dntp_conc=dntp_conc,
        dna_conc=dna_conc,
        tm_method='santalucia',
        salt_corrections_method='santalucia',
    )

@lru_cache(maxsize=None)
def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return str(Seq(seq).reverse_complement())


def gc_decimal(seq: str) -> float:
    """Return GC content of a sequence as a decimal."""
    s = seq.upper()
    return (s.count('G') + s.count('C')) / len(s)


def gc_percent(seq: str) -> float:
    """Return GC content of a sequence as a percentage (0–100)."""
    return round(gc_decimal(seq) * 100, 1)


def clean_input_sequence(seq: str) -> str:
    """
    Cleans an input sequence string, supporting raw ACGT, FASTA, and GenBank formats.
    """
    seq_stripped = seq.strip()
    
    # Check for GenBank format
    if "LOCUS" in seq_stripped and "ORIGIN" in seq_stripped:
        match = re.search(r'ORIGIN([\s\S]*?)(?://|$)', seq_stripped)
        if match:
            seq_stripped = match.group(1)
            
    # Check for FASTA format
    lines = seq_stripped.splitlines()
    if any(l.strip().startswith('>') for l in lines):
        seq_lines = []
        in_first = False
        extra_entries = 0
        for line in lines:
            stripped = line.strip()
            if stripped.startswith('>'):
                if in_first:
                    extra_entries += 1
                    break
                in_first = True
            elif in_first:
                seq_lines.append(stripped)
        if extra_entries > 0:
            _warnings.warn(
                f"Multi-entry FASTA detected ({extra_entries + 1} entries), Only the first entry will be used",
                UserWarning,
                stacklevel=2,
            )
        seq_stripped = ''.join(seq_lines)
        
    # Strip whitespace, digits (GenBank line numbers), and newlines
    return re.sub(r'[\s0-9]', '', seq_stripped).upper()


# SECTION 2: SEQUENCE QUALITY CHECKS

@lru_cache(maxsize=None)
def check_overlap_gc(overlap_seq: str,
                     gc_range: Tuple[float, float]) -> List[str]:
    """
    Check whether the overlap GC content is within the acceptable range.

    Args:
        overlap_seq : the overlap DNA sequence
        gc_range    : (min_gc%, max_gc%) acceptable range

    Returns:
        List of warning strings (empty = no issues)
    """
    gc = gc_percent(overlap_seq)
    min_gc, max_gc = gc_range
    warnings = []
    if gc > max_gc:
        warnings.append(f"Overlap GC% is {gc}% which is ABOVE the maximum {max_gc}%, risk of formation of secondary structures")
    elif gc < min_gc:
        warnings.append(f"Overlap GC% is {gc}% which is BELOW the minimum {min_gc}%, A-T rich overlaps may anneal unreliably during assembly")
    return warnings

@lru_cache(maxsize=None)
def check_homopolymer(seq: str, max_run: int) -> List[str]:
    """
    Detect homopolymer runs in a sequence.

    Args:
        seq     : DNA sequence to scan
        max_run : minimum run length to flag

    Returns:
        List of warning strings (empty = no issues)
    """
    warnings = []
    pattern = re.compile(r'(.)\1{' + str(max_run - 1) + r',}')
    for match in pattern.finditer(seq.upper()):
        run = match.group()
        warnings.append(f"Homopolymer run: '{run}' ({len(run)} bp of '{run[0]}') at position {match.start()}–{match.end()-1}, risk of polymerase slippage")
    return warnings

@lru_cache(maxsize=None)
def check_dinucleotide_repeat(seq: str, max_repeats: int) -> List[str]:
    """
    Detect consecutive dinucleotide repeats in a sequence.

    Args:
        seq         : DNA sequence to scan
        max_repeats : minimum number of consecutive repeats to flag

    Returns:
        List of warning strings (empty = no issues)
    """
    warnings = []
    pattern = re.compile(r'((.)(?!\2).)\1{' + str(max_repeats - 1) + r',}')
    for match in pattern.finditer(seq.upper()):
        repeat_unit = match.group(1)
        full_match  = match.group()
        n_repeats   = len(full_match) // 2
        warnings.append(f"Dinucleotide repeat: '{repeat_unit}' repeated {n_repeats}x, ('{full_match}') at position {match.start()}–{match.end()-1}, risk of polymerase slippage")
    return warnings

def check_repeat_vicinity(
    fragment: str,
    overlap_start: int,
    overlap_end: int,
    buffer: int,
    max_homopolymer: int,
    max_dinucleotide: int,
    overlap_hp_results: Optional[List[str]] = None,
    overlap_dn_results: Optional[List[str]] = None,
) -> List[str]:
    """
    Scan a buffered zone around the overlap for homopolymer / dinucleotide
    repeats that are adjacent to but not fully inside the overlap.

    Args:
        fragment           : the full fragment sequence
        overlap_start      : start position of the overlap
        overlap_end        : end position of the overlap
        buffer             : number of extra bases to scan on each side
        max_homopolymer    : homopolymer run threshold
        max_dinucleotide   : repeat count threshold
        overlap_hp_results : pre-computed homopolymer check on the overlap (optional)
        overlap_dn_results : pre-computed dinucleotide check on the overlap (optional)

    Returns:
        List of proximity warning strings (empty = no issues)
    """
    frag_len   = len(fragment)
    zone_start = max(0, overlap_start - buffer)
    zone_end   = min(frag_len, overlap_end + buffer)
    zone_seq   = fragment[zone_start: zone_end]

    zone_hp = check_homopolymer(zone_seq, max_homopolymer)
    zone_dn = check_dinucleotide_repeat(zone_seq, max_dinucleotide)

    # Use pre-computed overlap results to avoid redundant re-computation
    if overlap_hp_results is None:
        overlap_hp_results = check_homopolymer(
            fragment[overlap_start:overlap_end], max_homopolymer)
    if overlap_dn_results is None:
        overlap_dn_results = check_dinucleotide_repeat(
            fragment[overlap_start:overlap_end], max_dinucleotide)

    warnings = []
    if zone_hp and not overlap_hp_results:
        for w in zone_hp:
            warnings.append(f"Proximity warning: {w}")
    if zone_dn and not overlap_dn_results:
        for w in zone_dn:
            warnings.append(f"Proximity warning: {w}")
    return warnings

@lru_cache(maxsize=None)
def check_gc_clamp(
    oligo_seq: str,
    min_gc: int,
    clamp_length: int,
) -> List[str]:
    """
    Check whether the 3' end of an oligo has sufficient (min_gc) GC content (GC clamp).

    Args:
        oligo_seq    : full oligo sequence
        min_gc       : minimum number of G/C required in the clamp region
        clamp_length : number of terminal bases to examine

    Returns:
        List of warning strings (empty = no issues)
    """
    if not oligo_seq or len(oligo_seq) < clamp_length:
        return []

    clamp   = oligo_seq[-clamp_length:].upper()
    gc_count = sum(1 for b in clamp if b in ('G', 'C'))

    if gc_count < min_gc:
        return [f"Weak GC clamp: 3' terminal {clamp_length} bases are '{clamp}' ({gc_count}/{min_gc} G/C required)"]
    return []

def evaluate_dna_sequence(
    seq: str,
    max_homopolymer: int = 4,
    max_dinucleotide: int = 3
) -> dict:
    """
    Evaluates a DNA sequence for structural complexities such as homopolymers, 
    dinucleotide repeats, and GC content. Provides a detailed breakdown of 
    specific repeat sequences found.
    """
    from collections import Counter
    
    seq_upper = seq.upper()
    
    # 1. GC content
    gc_perc = gc_percent(seq_upper)
    
    # 2. Homopolymers
    hp_pattern = re.compile(r'(.)\1{' + str(max_homopolymer - 1) + r',}')
    hp_matches = [m.group() for m in hp_pattern.finditer(seq_upper)]
    hp_details = dict(Counter(hp_matches))
    
    # 3. Dinucleotide repeats
    dn_pattern = re.compile(r'((.)(?!\2).)\1{' + str(max_dinucleotide - 1) + r',}')
    dn_matches = [m.group() for m in dn_pattern.finditer(seq_upper)]
    dn_details = dict(Counter(dn_matches))
    
    return {
        "gc_percent": gc_perc,
        "homopolymers": {
            "total_count": len(hp_matches),
            "details": hp_details
        },
        "dinucleotide_repeats": {
            "total_count": len(dn_matches),
            "details": dn_details
        }
    }


# SECTION 2a: SLIDING WINDOW STRUCTURE CHECKS (HAIRPIN + HOMODIMER, ANY LENGTH)

@lru_cache(maxsize=None)
def _sliding_windows(seq: str) -> List[Tuple[int, str]]:
    """
    Generate 60-bp windows for hairpin/homodimer screening.

    For oligos > LONG_OLIGO_TERMINAL_ONLY_THRESHOLD (200 bp), only the first
    and last 60 bp are returned (the functional priming termini).  Internal
    secondary structures in long PCA oligos are transient and resolved during
    polymerase extension, and primer3's nearest-neighbor model cannot capture
    true long-range intramolecular folding anyway.

    For shorter oligos (61–200 bp), overlapping 60-bp windows with a 30-bp
    step are generated across the full length, plus a 3'-anchored window.

    Args:
        seq: full oligo sequence

    Returns:
        List of (start_position, 60-bp chunk) tuples
    """
    n   = len(seq)
    wnd = SLIDING_WINDOW_SIZE

    # Long oligos: terminal-only mode (first 60 bp + last 60 bp)
    if n > LONG_OLIGO_TERMINAL_ONLY_THRESHOLD:
        return [(0, seq[:wnd]), (n - wnd, seq[-wnd:])]

    # Standard sliding window for moderate-length oligos (61–200 bp)
    stp = SLIDING_WINDOW_STEP
    starts = list(range(0, n - wnd + 1, stp))
    last = n - wnd
    if not starts or starts[-1] != last:
        starts.append(last)
    return [(s, seq[s:s + wnd]) for s in starts]

@lru_cache(maxsize=None)
def _hairpin_check_sliding(
    seq: str,
    dg_threshold: float,
    mv_conc: float = DEFAULT_MV_CONC,
    dv_conc: float = DEFAULT_DV_CONC,
    dntp_conc: float = DEFAULT_DNTP_CONC,
    dna_conc: float = DEFAULT_DNA_CONC,
) -> Tuple[List[str], float]:
    """
    Checks an oligo of any length for hairpin formation using primer3 and flags it if dG <= dg_threshold.
    Uses 60-bp sliding windows for sequences > 60 bp, evaluating the worst-case window.

    Args:
        seq          : full oligo sequence
        dg_threshold : ΔG cutoff in kcal/mol
        mv_conc / dv_conc / dntp_conc / dna_conc : reaction conditions (mM/nM)

    Returns:
        (warnings: List[str], severity: float)
        severity = |dG - threshold| in kcal/mol; 0.0 if no hairpin detected
    """
    if not seq:
        return [], 0.0

    chunks = [(0, seq)] if len(seq) <= SLIDING_WINDOW_SIZE else _sliding_windows(seq)
    worst_dg   = 0.0
    worst_warn: List[str] = []
    worst_sev  = 0.0

    total = len(chunks)
    for idx, (start, chunk) in enumerate(chunks):
        try:
            r    = primer3.calc_hairpin(chunk, mv_conc=mv_conc, dv_conc=dv_conc,
                                        dntp_conc=dntp_conc, dna_conc=dna_conc)
            dg_k = r.dg / 1000.0 # convert result from cal/mol to kcal/mol
            tm_k = r.tm
        except Exception:
            continue
        if dg_k < worst_dg:
            worst_dg = dg_k
            if dg_k <= dg_threshold:
                loc = (f", window {idx+1}/{total} (bases {start}–{start+len(chunk)})" 
                       if len(seq) > SLIDING_WINDOW_SIZE else "")
                worst_sev  = abs(dg_k - dg_threshold)
                worst_warn = [f"Hairpin detected: dG = {dg_k:.2f} kcal/mol [Tm = {tm_k:.1f}°C] (threshold {dg_threshold:.1f}){loc}"]

    return worst_warn, worst_sev

@lru_cache(maxsize=None)
def _homodimer_check_sliding(
    seq: str,
    dg_threshold: float,
    mv_conc: float = DEFAULT_MV_CONC,
    dv_conc: float = DEFAULT_DV_CONC,
    dntp_conc: float = DEFAULT_DNTP_CONC,
    dna_conc: float = DEFAULT_DNA_CONC,
) -> Tuple[List[str], float]:
    """
    Checks an oligo of any length for homodimer formation using primer3 and flags it if dG <= dg_threshold.
    Uses 60-bp sliding windows for sequences > 60 bp, evaluating the worst-case window.

    Args:
        seq          : full oligo sequence
        dg_threshold : ΔG cutoff in kcal/mol
        mv_conc / dv_conc / dntp_conc / dna_conc : reaction conditions (mM/nM)

    Returns:
        (warnings: List[str], severity: float)
        severity = |dG - threshold| in kcal/mol; 0.0 if no homodimer detected
    """
    if not seq:
        return [], 0.0

    chunks = [(0, seq)] if len(seq) <= SLIDING_WINDOW_SIZE else _sliding_windows(seq)
    worst_dg   = 0.0
    worst_warn: List[str] = []
    worst_sev  = 0.0

    total = len(chunks)
    for i, (start, chunk) in enumerate(chunks):
        try:
            r    = primer3.calc_homodimer(chunk, mv_conc=mv_conc, dv_conc=dv_conc,
                                          dntp_conc=dntp_conc, dna_conc=dna_conc)
            dg_k = r.dg / 1000.0 # convert result from cal/mol to kcal/mol
            tm_k = r.tm
        except Exception:
            continue
        if dg_k < worst_dg:
            worst_dg = dg_k
            if dg_k <= dg_threshold:
                loc = (f", window {i+1}/{total} (bases {start}–{start+len(chunk)})"
                       if len(seq) > SLIDING_WINDOW_SIZE else "")
                worst_sev  = abs(dg_k - dg_threshold)
                worst_warn = [f"Homodimer risk: dG = {dg_k:.2f} kcal/mol [Tm = {tm_k:.1f}°C] (threshold {dg_threshold:.1f}){loc}"]

    return worst_warn, worst_sev

# SECTION 2b: MISPRIME RISK MAP & MISPRIME CHECK

def _build_tip_index(seq: str, probe_len: int, tip_length: int) -> dict:
    """
    Build a tip -> probe window start position lookup index for fast misprime candidate filtering.
    Maps each probe-window tip to all start positions of the probe window in seq.
    
    Args:
        seq        : sequence to index (top strand or reverse complement)
        probe_len  : probe window length
        tip_length : number of bp at the 3' end of each window used as the lookup key
    
    Returns:
        a dictionary of tip sequences mapped to the starting positions of its probe window
    """
    index = defaultdict(list)
    n = len(seq)
    for i in range(n - probe_len + 1):
        tip = seq[i + probe_len - tip_length: i + probe_len]
        index[tip].append(i)
    return dict(index)

def _build_misprime_risk_map(
    sequence: str,
    kmer_size: int = DEFAULT_MISPRIME_KMER_SIZE,
) -> Tuple[List[List[int]], List[List[int]]]:
    """
    Build per-position misprime risk maps for a sequence.
    
    Args:
        sequence  : DNA sequence
        kmer_size : 3' seed k-mer length for misprime detection
    
    Returns:
        (fwd_map, rev_map): Tuple[List[List[int]], List[List[int]]] of length n
        fwd_map[i]: annealing sites for a forward oligo with 3' end at position i
        rev_map[i]: annealing sites for a reverse oligo with top-strand start at position i
    """
    n        = len(sequence)
    full_seq = sequence.upper()

    # top-strand k-mer index: kmer → [start positions]
    top_idx: Dict[str, List[int]] = defaultdict(list)
    for i in range(n - kmer_size + 1):
        top_idx[full_seq[i: i + kmer_size]].append(i)
    top_idx = dict(top_idx)

    def _match_positions(probe_kmer: str) -> List[int]:
        """
        Return all top-strand positions where probe_kmer can anneal on both top and bottom strands
        """
        positions: List[int] = []

        # top-strand anneal
        for site in top_idx.get(probe_kmer, []):
            positions.append(site)

        # bottom-strand anneal
        rc_probe = reverse_complement(probe_kmer)
        for site in top_idx.get(rc_probe, []):
            positions.append(site)

        return positions

    # fwd_map: stores probe 12mer for a forward oligo and its top strand end positions
    fwd_map: List[List[int]] = [[] for _ in range(n)]
    for pos in range(kmer_size - 1, n):
        probe = full_seq[pos - kmer_size + 1: pos + 1]
        fwd_map[pos] = _match_positions(probe)

    # rev_map: stores probe 12mer for a reverse oligo and its top strand start positions
    rev_map: List[List[int]] = [[] for _ in range(n)]
    for pos in range(n - kmer_size + 1):
        probe = reverse_complement(full_seq[pos: pos + kmer_size])
        rev_map[pos] = _match_positions(probe)

    return fwd_map, rev_map

def check_misprime_all(
    records: List[dict],
    fragments: List[str],
    misprime_length: int,
    tip_length: int,
    max_mismatches: int,
    frag_global_starts: Optional[List[int]] = None,
) -> Tuple[List[bool], List[Optional[str]]]:
    """
    Post-tiling misprime check for each oligo's 3' probe against its own fragment sequence (both strands) excluding its own site
    
    Args:
        records             : list of oligo record dicts
        fragments           : list of fragment sequences
        misprime_length     : bp from the 3' end to use as the probe
        tip_length          : bp at the very 3' tip that must match exactly
        max_mismatches      : max mismatches allowed in the non-tip probe region
        frag_global_starts  : global start position of each fragment (accounts for junction overlaps).
                            If None, falls back to cumulative sum of fragment lengths.
    
    Returns:
        a list of misprime flags and misprime warnings (empty = no issues)
    """
    
    if frag_global_starts is None:
        _starts = []
        cumulative = 0
        for frag in fragments:
            _starts.append(cumulative)
            cumulative += len(frag)
        frag_global_starts = _starts

    # Build top + rc tip indexes from each fragment's own sequence
    frag_top_indexes = {}
    frag_rc_indexes  = {}
    for frag_idx, frag_seq in enumerate(fragments):
        frag_up = frag_seq.upper()
        frag_rc = reverse_complement(frag_up)
        frag_top_indexes[frag_idx] = _build_tip_index(frag_up, misprime_length, tip_length)
        frag_rc_indexes[frag_idx]  = _build_tip_index(frag_rc, misprime_length, tip_length)

    misprime_flags: List[bool]           = []
    misprime_warns: List[Optional[str]]  = []

    for row in records:
        oligo_seq  = str(row['oligo_sequence']).upper()
        frag_idx   = int(row['fragment']) - 1
        frag_start = int(row['frag_start'])
        frag_end   = int(row['frag_end'])

        # Report in global coordinates for clarity
        g_start = frag_global_starts[frag_idx]

        probe         = oligo_seq[-misprime_length:]
        probe_len     = len(probe)
        tip           = probe[-tip_length:]
        non_tip_probe = probe[:-tip_length]

        # Safe zone in fragment-local coordinates
        frag_seq   = fragments[frag_idx].upper()
        frag_rc    = reverse_complement(frag_seq)
        fseq_n     = len(frag_seq)
        safe_start  = max(0, frag_start - MISPRIME_SAFE_BUFFER)
        safe_end    = min(fseq_n, frag_end + MISPRIME_SAFE_BUFFER)

        hits = []
        top_index = frag_top_indexes[frag_idx]
        rc_index  = frag_rc_indexes[frag_idx]

        strand_configs = [
            ('top strand',    top_index, frag_seq, False),
            ('bottom strand', rc_index,  frag_rc,  True),
        ]

        for strand_label, index, search_seq, is_rc in strand_configs:
            candidate_positions = index.get(tip, [])
            for i in candidate_positions:
                if is_rc:
                    top_start = fseq_n - i - probe_len
                    top_end   = fseq_n - i
                else:
                    top_start = i
                    top_end   = i + probe_len

                # Skip the safe zone
                if not (top_end <= safe_start or top_start >= safe_end):
                    continue

                non_tip_window = search_seq[i: i + probe_len - tip_length]
                mismatches = sum(a != b for a, b in zip(non_tip_probe, non_tip_window))

                if mismatches <= max_mismatches:
                    hits.append((strand_label,
                                 g_start + top_start,
                                 g_start + top_end,
                                 mismatches))

        if hits:
            hit_descriptions = [
                f"{strand} pos {hs}–{he} ({mm} mismatch{'es' if mm != 1 else ''})"
                for strand, hs, he, mm in hits[:3]
            ]
            extra = f" (+{len(hits)-3} more)" if len(hits) > 3 else ""
            misprime_flags.append(True)
            misprime_warns.append(f"Misprime risk: 3' probe ({probe_len} bp) matches {len(hits)} unintended site(s), {', '.join(hit_descriptions)}{extra}")
        else:
            misprime_flags.append(False)
            misprime_warns.append(None)

    return misprime_flags, misprime_warns

# SECTION 2c: HETERODIMER RC K-MER RISK MAP & STRUCTURAL CHECK

def _build_hetdimer_risk_map(
    fragment: str,
    kmer_size: int = DEFAULT_HETDIMER_KMER_SIZE,
) -> Dict[str, List[int]]:
    """
    Build a top-strand k-mer position index for fragment-local heterodimer risk assessment.

    Args:
        fragment  : DNA fragment sequence (top strand)
        kmer_size : k-mer seed length for complementarity search

    Returns:
        a dictionary of k-mers and all of their start positions in the fragment
    """
    frag = fragment.upper()
    n    = len(frag)
    index: Dict[str, List[int]] = {}
    for i in range(n - kmer_size + 1):
        kmer = frag[i:i + kmer_size]
        if kmer in index:
            index[kmer].append(i)
        else:
            index[kmer] = [i]
    return index

def check_heterodimer(
    seq1: str,
    seq2: str,
    dg_threshold: float,
    mv_conc: float = DEFAULT_MV_CONC,
    dv_conc: float = DEFAULT_DV_CONC,
    dntp_conc: float = DEFAULT_DNTP_CONC,
    dna_conc: float = DEFAULT_DNA_CONC,
) -> Optional[str]:
    """
    Check whether two oligos form a heterodimer using primer3

    Args:
        seq1          : sequence of the first oligo
        seq2          : sequence of the second oligo
        dg_threshold  : ΔG cutoff in kcal/mol
        mv_conc / dv_conc / dntp_conc / dna_conc : reaction conditions (mM/nM)

    Returns:
        Warning string if heterodimer is stable enough to flag, else None.
    """
    if not seq1 or not seq2:
        return None

    s1 = seq1[-MAX_PRIMER3_SEQ_LENGTH:] if len(seq1) > MAX_PRIMER3_SEQ_LENGTH else seq1
    s2 = seq2[-MAX_PRIMER3_SEQ_LENGTH:] if len(seq2) > MAX_PRIMER3_SEQ_LENGTH else seq2

    result = primer3.calc_heterodimer(s1.upper(), s2.upper(), mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=dna_conc)
    dg = result.dg / 1000.0 # convert result from cal/mol to kcal/mol
    tm = result.tm

    if dg <= dg_threshold:
        trunc_note = (" [3'-terminal 60 bp only]" if (len(seq1) > MAX_PRIMER3_SEQ_LENGTH or len(seq2) > MAX_PRIMER3_SEQ_LENGTH) else "")
        return (f"Heterodimer risk: dG = {dg:.2f} kcal/mol [Tm = {tm:.1f}°C] (threshold {dg_threshold:.1f}){trunc_note}")
    return None

def check_heterodimers_all(
    records: List[dict],
    dg_threshold: float = DEFAULT_HETERODIMER_DG_THRESHOLD,
    neighbor_range: int  = DEFAULT_HETERODIMER_NEIGHBOR_RANGE,
    dg_buffer: float = DEFAULT_HETERODIMER_DG_BUFFER,
    mv_conc: float = DEFAULT_MV_CONC,
    dv_conc: float = DEFAULT_DV_CONC,
    dntp_conc: float = DEFAULT_DNTP_CONC,
    dna_conc: float = DEFAULT_DNA_CONC,
) -> Tuple[List[bool], List[Optional[str]]]:
    """
    Post-tiling Fragment-local heterodimer check for each oligo against ±neighbor_range non-adjacent neighbors (skipping ±1).

    Args:
        records        : list of all oligo record dicts
        dg_threshold   : ΔG cutoff for flagging (kcal/mol)
        neighbor_range : check up to ±N neighbors, skipping ±1
        dg_buffer      : buffer applied to threshold before flagging
        mv_conc / dv_conc / dntp_conc / dna_conc : reaction conditions (mM/nM)

    Returns:
        hetdimer_flags: list[bool], hetdimer_warns: list[str | None]
    """
    n        = len(records)
    seqs     = [str(r['oligo_sequence']).upper() for r in records]
    frag_ids = [r.get('fragment', 0) for r in records] # to skip cross-fragment checks

    hetdimer_flags: List[bool]          = []
    hetdimer_warns: List[Optional[str]] = []

    for i in range(n):
        hits = []

        for offset in range(2, neighbor_range + 1):
            for j in [i - offset, i + offset]:
                if j < 0 or j >= n:
                    continue   # outside of list bounds
                if frag_ids[i] != frag_ids[j]:
                    continue   # cross-fragment pair
                warn = check_heterodimer(seqs[i], seqs[j], dg_threshold - dg_buffer, mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=dna_conc)
                if warn:
                    neighbor_name = records[j].get('oligo_name', f'oligo_{j+1}')
                    hits.append(
                        f"vs {neighbor_name} (offset {'+' if j>i else ''}{j-i}): {warn}"
                    )

        if hits:
            hetdimer_flags.append(True)
            hetdimer_warns.append('; '.join(hits))
        else:
            hetdimer_flags.append(False)
            hetdimer_warns.append(None)

    return hetdimer_flags, hetdimer_warns

# SECTION 3: FRAGMENTATION

def _score_cut_point(sequence: str, pos: int, window: int,
                     max_homopolymer: int, max_dinucleotide: int,
                     repeat_penalty: float = 1000.0) -> float:
    """
    Score a candidate fragment cut point, lower = better

    Args:
        sequence         : full DNA sequence
        pos              : candidate cut position
        window           : bp window around cut to examine
        max_homopolymer  : homopolymer run length threshold
        max_dinucleotide : dinucleotide repeat length threshold
        repeat_penalty   : penalty per repeats in the window

    Returns:
        Penalty Score (float)
    """
    half  = window // 2
    n     = len(sequence)
    start = max(0, pos - half)
    end   = min(n, pos + half)
    region = sequence[start:end]

    gc_score   = abs(gc_decimal(region) * 100 - 50.0)
    hp_penalty = len(check_homopolymer(region, max_homopolymer)) * repeat_penalty
    dn_penalty = len(check_dinucleotide_repeat(region, max_dinucleotide)) * repeat_penalty

    return gc_score + hp_penalty + dn_penalty

def _fragment_sequence(
    sequence: str,
    frag_range: Tuple[int, int],
    junction_overlap: int = 0,
    cut_search_window: int = DEFAULT_CUT_SEARCH_WINDOW,
    max_homopolymer: int   = DEFAULT_MAX_HOMOPOLYMER_RUN,
    max_dinucleotide: int  = DEFAULT_MAX_DINUCLEOTIDE_REPEAT,
    cut_score_window: int  = 20,
    cut_score_repeat_penalty: float = 1000.0,
) -> List[str]:
    """
    Divide the full sequence into fragments within frag_range optimizing for to minimize extreme GC% and repeats
    
    Args:
        sequence                 : full DNA sequence
        frag_range               : fragment length bounds
        junction_overlap         : bp shared at each fragment boundary for assembly
        cut_search_window        : bp to search on either side of the cut point
        max_homopolymer          : homopolymer run length threshold
        max_dinucleotide         : dinucleotide repeat length threshold
        cut_score_window         : bp window used by _score_cut_point
        cut_score_repeat_penalty : penalty per repeat element at the cut site

    Returns:
        A list of fragments with overlapping junctions
    """
    min_f, max_f = frag_range
    target = (min_f + max_f) // 2
    n      = len(sequence)

    if n <= max_f:
        return [sequence]

    # find lowest penalty cut points
    cut_points = []
    pos = 0
    while pos < n:
        remaining = n - pos
        if remaining <= max_f:
            cut_points.append(pos)
            break
        cut_points.append(pos)

        naive_cut = pos + target
        search_lo = max(pos + min_f, naive_cut - cut_search_window)
        search_hi = min(pos + max_f, naive_cut + cut_search_window, n)

        best_cut   = naive_cut
        best_score = _score_cut_point(sequence, naive_cut, cut_score_window,
                                      max_homopolymer, max_dinucleotide,
                                      cut_score_repeat_penalty)
        for candidate in range(search_lo, search_hi + 1):
            score = _score_cut_point(sequence, candidate, cut_score_window,
                                     max_homopolymer, max_dinucleotide,
                                     cut_score_repeat_penalty)
            if score < best_score:
                best_score = score
                best_cut   = candidate

        pos = best_cut

    # Merge tiny terminal fragment into the previous one
    if len(cut_points) >= 2:
        last_start    = cut_points[-1]
        terminal_len  = n - last_start
        if terminal_len < min_f:
            _warnings.warn(
                f"Terminal fragment is {terminal_len} bp (< min {min_f} bp), merging into previous fragment to form a {n - cut_points[-2]} bp fragment (max={max_f})",
                UserWarning,
                stacklevel=3,
            )
            cut_points.pop()

    # add Junction Overlaps to each fragment (except last one)
    fragments = []
    for idx, start in enumerate(cut_points):
        if idx < len(cut_points) - 1:
            raw_end = cut_points[idx + 1]
            end     = min(n, raw_end + junction_overlap)
        else:
            end = n
        fragments.append(sequence[start:end])

    return fragments

def _fragment_sequence_golden_gate(
    sequence: str,
    frag_range: Tuple[int, int],
    cut_search_window: int,
    max_homopolymer: int,
    max_dinucleotide: int,
    cut_score_window: int,
    cut_score_repeat_penalty: float,
    enzyme_name: str,
) -> Tuple[List[str], List[str]]:
    min_f, max_f = frag_range
    target = (min_f + max_f) // 2
    n = len(sequence)

    if n <= max_f:
        return [sequence], []

    def _rc(seq: str) -> str:
        return reverse_complement(seq)

    def _hamming(a: str, b: str) -> int:
        return sum(x != y for x, y in zip(a, b))

    def _is_valid_overhang(ov: str, chosen: List[str]) -> bool:
        if len(ov) != 4: return False
        if _rc(ov) == ov: return False # palindrome
        for homop in ['AAAA', 'TTTT', 'CCCC', 'GGGG']:
            if ov == homop: return False

        rec = GOLDEN_GATE_ENZYMES[enzyme_name]["recognition"]
        if ov.startswith(rec[:4]) or _rc(ov).startswith(rec[:4]):
            return False

        # orthogonality: differ from every chosen overhang, and its reverse complement,
        # by at least GOLDEN_GATE_MIN_OVERHANG_DISTANCE bases
        for prev in chosen:
            if (_hamming(ov, prev) < GOLDEN_GATE_MIN_OVERHANG_DISTANCE
                    or _hamming(ov, _rc(prev)) < GOLDEN_GATE_MIN_OVERHANG_DISTANCE):
                return False

        return True

    cut_points = []
    chosen_overhangs = []
    pos = 0
    while pos < n:
        remaining = n - pos
        if remaining <= max_f:
            cut_points.append(pos)
            break
        cut_points.append(pos)

        naive_cut = pos + target
        search_lo = max(pos + min_f, naive_cut - cut_search_window)
        search_hi = min(pos + max_f, naive_cut + cut_search_window, n)

        best_cut = None
        best_score = float('inf')
        
        for candidate in range(search_lo, search_hi + 1):
            if candidate + 4 > n:
                continue
            ov = sequence[candidate:candidate+4]
            if not _is_valid_overhang(ov, chosen_overhangs):
                continue
                
            score = _score_cut_point(sequence, candidate, cut_score_window,
                                     max_homopolymer, max_dinucleotide,
                                     cut_score_repeat_penalty)
            if score < best_score:
                best_score = score
                best_cut = candidate
                
        if best_cut is None:
            raise ValueError(
                f"Could not find a valid native Golden Gate overhang in the window [{search_lo}, {search_hi}] "
                f"after {len(chosen_overhangs)} overhang(s) were chosen: every candidate 4-mer is a palindrome, a homopolymer, "
                f"starts like the enzyme recognition site, or differs by fewer than {GOLDEN_GATE_MIN_OVERHANG_DISTANCE} bases "
                f"from an earlier overhang or its reverse complement. Use larger fragments (fewer junctions) or Gibson Assembly instead."
            )
            
        pos = best_cut
        chosen_overhangs.append(sequence[best_cut:best_cut+4])

    if len(cut_points) >= 2:
        last_start = cut_points[-1]
        terminal_len = n - last_start
        if terminal_len < min_f:
            _warnings.warn(
                f"Terminal fragment is {terminal_len} bp (< min {min_f} bp), merging into previous fragment to form a {n - cut_points[-2]} bp fragment (max={max_f})",
                UserWarning,
                stacklevel=3,
            )
            cut_points.pop()
            chosen_overhangs.pop()

    fragments = []
    for idx, start in enumerate(cut_points):
        if idx < len(cut_points) - 1:
            end = cut_points[idx + 1] + 4
        else:
            end = n
        fragments.append(sequence[start:end])

    return fragments, chosen_overhangs

# SECTION 4a: RECORD BUILDER

def _make_record(
    frag_idx: int, oligo_idx: int, direction: str,
    frag_start: int, frag_end: int,
    oligo_seq: str, overlap_seq: Optional[str],
    overlap_tm: Optional[float], warnings: List[str], is_last: bool,
    tm_flag: bool = False,
    gc_flag: bool = False,
    homopolymer_flag: bool = False,
    dinucleotide_flag: bool = False,
    hairpin_flag: bool = False,
    homodimer_flag: bool = False,
    gc_clamp_flag: bool = False,
    heterodimer_flag: bool = False,
    body_homopolymer_flag: bool = False,
    body_dinucleotide_flag: bool = False,
) -> dict:
    """
    Build a single oligo record
    """
    dir_tag = 'F' if direction == 'forward' else 'R'
    return {
        # Oligo information
        'oligo_name':               f"Frag{frag_idx+1:02d}_Oligo{oligo_idx+1:02d}_{dir_tag}",
        'fragment':                 frag_idx + 1,
        'oligo_in_frag':            oligo_idx + 1,
        'direction':                direction,

        # Position on fragment
        'frag_start':               frag_start,
        'frag_end':                 frag_end,

        # Oligo properties
        'oligo_length':             frag_end - frag_start,
        'oligo_gc_%':               gc_percent(oligo_seq),
        'oligo_sequence':           oligo_seq,

        # Overlap properties
        'overlap_sequence':         overlap_seq if overlap_seq else 'N/A',
        'overlap_length':           len(overlap_seq) if overlap_seq else None,
        'overlap_tm_C':             overlap_tm if not is_last else None,
        'overlap_gc_%':             gc_percent(overlap_seq) if overlap_seq else None,

        # Quality flags
        'flag_tm_out_of_range':     tm_flag,
        'flag_overlap_gc':          gc_flag,
        'flag_homopolymer':         homopolymer_flag,
        'flag_dinucleotide_repeat': dinucleotide_flag,
        'flag_hairpin':             hairpin_flag,
        'flag_homodimer':           homodimer_flag,
        'flag_gc_clamp':            gc_clamp_flag,
        'flag_heterodimer':         heterodimer_flag,
        'flag_body_homopolymer':    body_homopolymer_flag,
        'flag_body_dinucleotide':   body_dinucleotide_flag,
        'flag_misprime':            False,

        # A* DP penalty (tiling-time cost only; post-tiling checks not included)
        'penalty_score':            0.0,

        # Warnings
        'is_last_in_frag':          is_last,
        'warnings':                 '; '.join(warnings) if warnings else 'OK',
    }

# SECTION 4b: GOLDEN GATE HELPER FUNCTION

def _scan_restriction_sites(sequence: str, enzyme_name: str) -> List[int]:
    """
    Return all positions where the enzyme recognition site (or its RC) appears in sequence

    Args:
        sequence    : full DNA sequence
        enzyme_name : key into GOLDEN_GATE_ENZYMES dict

    Returns:
        a list of sorted 0-based start positions for an enzyme's recognition site (or its RC), Returns [] if no sites found
    """
    rec    = GOLDEN_GATE_ENZYMES[enzyme_name]["recognition"]
    rc_rec = GOLDEN_GATE_ENZYMES[enzyme_name]["rc_recognition"]
    seq    = sequence.upper()
    positions = []
    for i in range(len(seq) - len(rec) + 1):
        if seq[i:i+len(rec)] in (rec, rc_rec):
            positions.append(i)
    return sorted(positions)

def _design_goldengate_overhangs(n_overhangs: int, enzyme_name: str) -> List[str]:
    """
    Design mutually orthogonal 4bp overhangs for Golden Gate Assembly.

    Args:
        n_overhangs : overhangs needed (= n_fragments + 1)
        enzyme_name : key into GOLDEN_GATE_ENZYMES dict
    
    Returns:
        a list of orthogonal 4-mer sequences overhangs
    
    Raises:
        ValueError if the orthogonal pool is exhausted before finding enough overhangs
    """
    from itertools import product as _product

    def _hamming(a: str, b: str) -> int:
        return sum(x != y for x, y in zip(a, b))

    def _rc(seq: str) -> str:
        return reverse_complement(seq)

    bases = 'ACGT'
    # generate all 4mer overhangs, deduplicate via RC pairs
    candidates = []
    seen = set()
    for combo in _product(bases, repeat=4):
        mer = ''.join(combo)
        if mer in seen:
            continue
        rc_mer = _rc(mer)
        if rc_mer == mer:
            # remove palindromes
            seen.add(mer)
            continue
        # keep lexicographically smaller
        keep = min(mer, rc_mer)
        if keep not in seen:
            candidates.append(keep)
            seen.add(keep)
            seen.add(_rc(keep))

    # remove homopolymers
    homopolymers = {'AAAA', 'TTTT', 'CCCC', 'GGGG'}
    candidates = [c for c in candidates if c not in homopolymers]

    # remove any overhang starting with the enzyme recognition sequence
    rec = GOLDEN_GATE_ENZYMES[enzyme_name]["recognition"]
    candidates = [
        c for c in candidates
        if not c.startswith(rec[:4]) and not _rc(c).startswith(rec[:4])
    ]

    # sorting (for reproducibility)
    candidates.sort()

    # greedy orthogonal selection
    selected: List[str] = []
    for cand in candidates:
        rc_cand = _rc(cand)
        ok = True
        for sel in selected:
            rc_sel = _rc(sel)
            # Hamming distance between cand and sel
            if _hamming(cand, sel) < GOLDEN_GATE_MIN_OVERHANG_DISTANCE:
                ok = False; break
            # cand vs RC(sel)
            if _hamming(cand, rc_sel) < GOLDEN_GATE_MIN_OVERHANG_DISTANCE:
                ok = False; break
            # RC(cand) vs sel
            if _hamming(rc_cand, sel) < GOLDEN_GATE_MIN_OVERHANG_DISTANCE:
                ok = False; break
        if ok:
            selected.append(cand)
            if len(selected) == n_overhangs:
                break

    if len(selected) < n_overhangs:
        raise ValueError(f"Golden Gate overhang pool exhausted: only {len(selected)} of {n_overhangs} orthogonal 4-bp overhangs found. Reduce fragment count (max: {MAX_GOLDENGATE_FRAGMENTS}).")
    return selected



def _design_fragment_table(
    fragments: List[str],
    junction_overlap: int,
    junction_tm_range: Tuple[float, float],
    tm_warning_buffer: float,
    gc_range: Tuple[float, float],
    max_homopolymer: int,
    max_dinucleotide: int,
    hairpin_dg_threshold: float,
    mv_conc: float,
    dv_conc: float,
    dntp_conc: float,
    dna_conc: float,
    frag_global_starts: Optional[List[int]] = None,
    assembly_method: str = DEFAULT_ASSEMBLY_METHOD,
) -> pd.DataFrame:
    """
    Build the fragment summary table: one row per fragment with metadata, junction sequence, and quality flags.
    
    Args:
        fragments          : fragment sequences 
        junction_overlap   : bp shared between consecutive fragments
        junction_tm_range  : acceptable Tm range for the junction overlap sequence
        tm_warning_buffer  : tolerance added to junction_tm_range bounds before flagging
        gc_range           : GC% bounds for junction quality check
        max_homopolymer    : homopolymer run length threshold
        max_dinucleotide   : dinucleotide repeat threshold
        hairpin_dg_threshold : ΔG threshold for hairpin detection
        mv_conc, dv_conc, dntp_conc, dna_conc : reaction conditions
        frag_global_starts : pre-computed global start per fragment; computed internally if None
        assembly_method    : 'stitching_pcr' | 'gibson' | 'golden_gate' (affects junction column names)
    
    Returns:
        pd.DataFrame with one row per fragment
    """
    if not fragments:
        return pd.DataFrame()

    # Compute global starts if not supplied
    if frag_global_starts is None:
        starts = []
        cumulative = 0
        for i, frag in enumerate(fragments):
            starts.append(cumulative)
            if i < len(fragments) - 1:
                cumulative += len(frag) - junction_overlap
            else:
                cumulative += len(frag)
        frag_global_starts = starts

    records = []

    for i, frag in enumerate(fragments):
        frag_num    = i + 1
        frag_len    = len(frag)
        g_start     = frag_global_starts[i]

        # true end of this fragment's unique region (excludes junction overlap)
        if i < len(fragments) - 1:
            g_end_unique = g_start + frag_len - junction_overlap
        else:
            g_end_unique = g_start + frag_len

        row: dict = {
            'fragment_num'       : frag_num,
            'fragment_length_bp' : frag_len,
            'global_start'       : g_start,
            'global_end'         : g_end_unique,
            'fragment_sequence'  : frag,
        }

        # Junction to next fragment
        if i < len(fragments) - 1 and junction_overlap > 0:
            frag_next = fragments[i + 1]

            junc_seq   = frag[-junction_overlap:]
            junc_start = g_start + frag_len - junction_overlap
            junc_end   = junc_start + len(junc_seq)

            # just a sanity check, junction must match the start of the next fragment
            expected_in_next    = frag_next[:junction_overlap]
            sequence_consistent = (junc_seq.upper() == expected_in_next.upper())

            tm      = calc_tm(junc_seq, mv_conc, dv_conc, dntp_conc, dna_conc)
            tm_ok   = (junction_tm_range[0] - tm_warning_buffer) <= tm <= (junction_tm_range[1] + tm_warning_buffer)
            gc_warn = check_overlap_gc(junc_seq, gc_range)
            hp_warn = check_homopolymer(junc_seq, max_homopolymer)
            dn_warn = check_dinucleotide_repeat(junc_seq, max_dinucleotide)

            fwd_primer     = junc_seq
            rev_primer     = reverse_complement(junc_seq)

            hp_fwd         = _hairpin_check_sliding(
                fwd_primer, hairpin_dg_threshold,
                mv_conc, dv_conc, dntp_conc, dna_conc,
            )[0] 
            hp_rev         = _hairpin_check_sliding(
                rev_primer, hairpin_dg_threshold,
                mv_conc, dv_conc, dntp_conc, dna_conc,
            )[0]
            hp_primer_warn = hp_fwd + hp_rev

            all_warnings = []
            if not tm_ok:
                direction_word = ('BELOW minimum' if tm < (junction_tm_range[0] - tm_warning_buffer)
                                  else 'ABOVE maximum')
                bound = (junction_tm_range[0] if tm < (junction_tm_range[0] - tm_warning_buffer)
                         else junction_tm_range[1])
                # Label the warning with the active assembly method context
                _method_label = (
                    "Gibson Assembly annealing (50°C isothermal)"
                    if assembly_method == "gibson"
                    else "stitching PCR"
                )
                all_warnings.append(
                    f"Junction Tm {tm:.1f}C is {direction_word} {bound}C, {_method_label} may not work at your intended annealing temperature"
                )
            all_warnings += gc_warn + hp_warn + dn_warn + hp_primer_warn
            if not sequence_consistent:
                all_warnings.append(
                    "WARNING: junction sequence does not match between this fragment "
                    "and the next, possible fragmentation bug"
                )

            # column names depend on assembly method
            if assembly_method == "gibson":
                _fwd_col = "fwd_gibson_overlap"
                _rev_col = "rev_gibson_overlap"
            elif assembly_method == "golden_gate":
                # Golden Gate: full primers are built separately with BsaI cassettes, these columns report the raw junction boundary sequences
                _fwd_col = "fwd_junction_boundary"
                _rev_col = "rev_junction_boundary"
            else:
                _fwd_col = "fwd_stitch_primer"
                _rev_col = "rev_stitch_primer"

            row.update({
                'junction_to_frag'      : frag_num + 1,
                'junction_start'        : junc_start,
                'junction_end'          : junc_end,
                'junction_length'       : len(junc_seq),
                'junction_sequence'     : junc_seq,
                _fwd_col                : fwd_primer,
                _rev_col                : rev_primer,
                'junction_tm_C'         : round(tm, 2),
                'junction_gc_%'         : gc_percent(junc_seq),
                'flag_tm_out_of_range'  : not tm_ok,
                'flag_overlap_gc'       : bool(gc_warn),
                'flag_homopolymer'      : bool(hp_warn),
                'flag_dinucleotide'     : bool(dn_warn),
                'flag_hairpin'          : bool(hp_primer_warn),
                'warnings'              : '; '.join(all_warnings) if all_warnings else 'OK',
            })
        else:
            # Last fragment (no junction)
            row.update({
                'junction_to_frag'      : None,
                'junction_start'        : None,
                'junction_end'          : None,
                'junction_length'       : None,
                'junction_sequence'     : None,
                'fwd_stitch_primer'     : None,
                'rev_stitch_primer'     : None,
                'junction_tm_C'         : None,
                'junction_gc_%'         : None,
                'flag_tm_out_of_range'  : False,
                'flag_overlap_gc'       : False,
                'flag_homopolymer'      : False,
                'flag_dinucleotide'     : False,
                'flag_hairpin'          : False,
                'warnings'              : 'OK',
            })

        records.append(row)

    # Cross-junction k-mer check for Gibson/stitching_pcr; for Golden Gate, overhang orthogonality (minimum Hamming distance) is enforced when the cut sites are chosen in _fragment_sequence_golden_gate()

    if assembly_method in ("gibson", "stitching_pcr") and len(records) > 1:
        _junc_data = [
            (r['junction_to_frag'], r['junction_sequence'])
            for r in records
            if r.get('junction_sequence') is not None
        ]
        _klen = MIN_RELIABLE_TM_LENGTH  # k-mer length = MIN_RELIABLE_TM_LENGTH
        _junc_kmers = [
            (frag_id, seq, {seq[k:k+_klen] for k in range(len(seq) - _klen + 1)})
            for frag_id, seq in _junc_data
        ]
        for _a in range(len(_junc_kmers)):
            for _b in range(_a + 1, len(_junc_kmers)):
                _id_a, _seq_a, _km_a = _junc_kmers[_a]
                _id_b, _seq_b, _km_b = _junc_kmers[_b]
                _shared = _km_a & _km_b
                if _shared:
                    _example = next(iter(_shared))
                    _warnings.warn(
                        f"Cross-junction similarity: junction {_id_a-1}→{_id_a} and {_id_b-1}→{_id_b} "
                        f"share {_klen}-mer {_example!r}. Risk of cross-hybridisation during assembly.",
                        UserWarning,
                        stacklevel=3,
                    )
                    # Also flag in the DataFrame row's warnings column
                    for r in records:
                        if r.get('junction_to_frag') in (_id_a, _id_b):
                            _prev = r.get('warnings', 'OK')
                            _note = (
                                f"Cross-junction k-mer match with junction "
                                f"{'%d→%d' % (_id_b-1, _id_b) if r['junction_to_frag'] == _id_a else '%d→%d' % (_id_a-1, _id_a)}: "
                                f"shared {_klen}-mer {_example!r}"
                            )
                            r['warnings'] = (_prev + '; ' + _note) if _prev != 'OK' else _note

    return pd.DataFrame(records)

# SECTION 5: OUTER PCR PRIMER DESIGN

def _design_outer_primers(
    fragments: List[str],
    primer_length_range: Tuple[int, int],
    primer_tm_range: Tuple[float, float],
    tm_warning_buffer: float,
    gc_range: Tuple[float, float],
    max_homopolymer: int,
    hairpin_dg_threshold: float,
    mv_conc: float,
    dv_conc: float,
    dntp_conc: float,
    dna_conc: float,
) -> pd.DataFrame:
    """
    Design outer PCR amplification primers for each fragment.

    Args:
        fragments           : list of fragment sequences to design primers for
        primer_length_range : (min_bp, max_bp)
        primer_tm_range     : (min_C, max_C) acceptable Tm range
        tm_warning_buffer   : tolerance outside primer_tm_range before flagging
        gc_range            : (min_%, max_%) acceptable GC content
        max_homopolymer     : homopolymer run threshold (flags 3' tail)
        hairpin_dg_threshold: ΔG threshold in kcal/mol
        mv_conc, dv_conc, dntp_conc, dna_conc : thermodynamic conditions

    Returns:
        pd.DataFrame with two rows per fragment (forward + reverse)
    """
    min_len, max_len = primer_length_range
    min_tm,  max_tm  = primer_tm_range
    min_gc,  max_gc  = gc_range
    opt_tm           = (min_tm + max_tm) / 2.0

    records = []

    for frag_idx, fragment in enumerate(fragments):
        frag_num = frag_idx + 1
        frag_len = len(fragment)

        for direction in ('forward', 'reverse'):
            dir_tag     = 'F' if direction == 'forward' else 'R'
            primer_name = f"Frag{frag_num:02d}_Outer_{dir_tag}"

            # Scan from the fragment edge inward; pick length with Tm closest to opt_tm.
            best_seq   = None
            best_tm_v  = None
            best_score = float('inf')

            for length in range(min_len, min(max_len, frag_len) + 1):
                cand  = (fragment[:length] if direction == 'forward'
                         else reverse_complement(fragment[-length:]))
                tm    = calc_tm(cand, mv_conc, dv_conc, dntp_conc, dna_conc)
                score = abs(tm - opt_tm)
                if score < best_score:
                    best_score = score
                    best_seq   = cand
                    best_tm_v  = tm

            if best_seq is None:
                records.append({
                    'primer_name'          : primer_name,
                    'fragment'             : frag_num,
                    'direction'            : direction,
                    'primer_sequence'      : 'N/A',
                    'primer_length'        : None,
                    'primer_tm_C'          : None,
                    'primer_gc_%'          : None,
                    'primer_penalty'       : None,
                    'flag_tm_out_of_range' : True,
                    'flag_gc'              : True,
                    'flag_homopolymer_3p'  : False,
                    'flag_hairpin'         : False,
                    'warnings'             : (
                        f"Fragment too short ({frag_len} bp) for primer length "
                        f"range [{min_len}, {max_len}] bp."
                    ),
                })
                continue

            primer_seq = best_seq
            primer_tm  = best_tm_v
            primer_gc  = gc_percent(primer_seq)

            hp_warn    = _hairpin_check_sliding(
                primer_seq, hairpin_dg_threshold,
                mv_conc, dv_conc, dntp_conc, dna_conc,
            )[0]
            tm_ok      = (min_tm - tm_warning_buffer) <= primer_tm <= (max_tm + tm_warning_buffer)
            gc_ok      = min_gc <= primer_gc <= max_gc
            tail       = primer_seq[-6:]
            hp_tail_ok = not check_homopolymer(tail, max_homopolymer)

            warns = []
            if not tm_ok:
                bound          = min_tm if primer_tm < (min_tm - tm_warning_buffer) else max_tm
                direction_word = ('BELOW minimum' if primer_tm < (min_tm - tm_warning_buffer)
                                  else 'ABOVE maximum')
                warns.append(
                    f"Primer Tm {primer_tm:.1f}\u00b0C is {direction_word} {bound}\u00b0C, "
                    f"PCR may not work at your intended annealing temperature"
                )
            if not gc_ok:
                warns.append(
                    f"Primer GC% {primer_gc:.1f}% is outside target range [{min_gc}%, {max_gc}%]"
                )
            if not hp_tail_ok:
                warns.append(
                    f"Homopolymer run at 3' end ('{tail}'), risk of polymerase slippage during PCR extension"
                )
            warns += hp_warn

            records.append({
                'primer_name'          : primer_name,
                'fragment'             : frag_num,
                'direction'            : direction,
                'primer_sequence'      : primer_seq,
                'primer_length'        : len(primer_seq),
                'primer_tm_C'          : round(primer_tm, 2),
                'primer_gc_%'          : round(primer_gc, 1),
                'primer_penalty'       : None,
                'flag_tm_out_of_range' : not tm_ok,
                'flag_gc'              : not gc_ok,
                'flag_homopolymer_3p'  : not hp_tail_ok,
                'flag_hairpin'         : bool(hp_warn),
                'warnings'             : '; '.join(warns) if warns else 'OK',
            })

    return pd.DataFrame(records)

# SECTION 6: TILING, HEAP-BASED A* SEARCH / DIJKSTRA ALGORITHM

def _scaled(base: float, severity: float, multiplier: float) -> float:
    """
    Scale a base DP penalty by how severe the quality violation is.
    Formula: base * (1 + multiplier * severity)
    """
    return base * (1.0 + multiplier * max(0.0, severity))

def _evaluate_placement(
    fragment: str, pos: int, oligo_len: int, ov_len: int, direction: str,
    tm_range: Tuple[float, float],
    tm_warning_buffer: float,
    gc_range: Tuple[float, float],
    max_homopolymer: int,
    max_dinucleotide: int,
    vicinity_buffer: int,
    hairpin_dg_threshold: float,
    mv_conc: float, dv_conc: float, dntp_conc: float, dna_conc: float,
    risk_map: Optional[Tuple[List[List[int]], List[List[int]]]] = None,
    frag_global_start: int = 0,
    penalties: Optional[dict] = None,
    homodimer_dg_threshold: float = DEFAULT_HOMODIMER_DG_THRESHOLD,
    gc_clamp_length: int = DEFAULT_GC_CLAMP_LENGTH,
    gc_clamp_min_gc: int = DEFAULT_GC_CLAMP_MIN_GC,
    misprime_kmer_size: int = DEFAULT_MISPRIME_KMER_SIZE,
    hetdimer_risk_map: Optional[Dict[str, List[int]]] = None,
    hetdimer_kmer_size: int = DEFAULT_HETDIMER_KMER_SIZE,
    gc_warning_buffer: float = DEFAULT_GC_WARNING_BUFFER,
    hairpin_dg_buffer: float = DEFAULT_HAIRPIN_DG_BUFFER,
    homodimer_dg_buffer: float = DEFAULT_HOMODIMER_DG_BUFFER,
) -> Tuple[float, dict]:
    """
    Compute the DP penalty and quality data for a single oligo placement

    Args:
        fragment              : fragment sequence being tiled
        pos                   : fragment-local oligo start (0-based)
        oligo_len             : oligo length in bp
        ov_len                : overlap length at the 3' end
        direction             : 'forward' or 'reverse'
        tm_range              : acceptable overlap Tm
        tm_warning_buffer     : tolerance outside tm_range before flagging
        gc_range              : acceptable overlap GC%
        max_homopolymer       : homopolymer run length threshold
        max_dinucleotide      : dinucleotide repeat threshold
        vicinity_buffer       : bp around overlap boundaries to check for repeats
        hairpin_dg_threshold  : ΔG threshold for hairpin flagging (kcal/mol)
        mv_conc, dv_conc, dntp_conc, dna_conc : reaction conditions
        risk_map              : (fwd_map, rev_map) from _build_misprime_risk_map(); None disables misprime penalty
        frag_global_start     : retained for API compatibility; unused (risk map is fragment-local)
        penalties             : penalty weight overrides; falls back to module-level defaults if None
        homodimer_dg_threshold: ΔG threshold for homodimer flagging (kcal/mol)
        gc_clamp_length       : number of 3' bases to check for GC clamp
        gc_clamp_min_gc       : minimum G/C bases required in gc_clamp_length
        misprime_kmer_size    : k-mer size for misprime risk lookup
        hetdimer_risk_map     : RC k-mer position index for heterodimer risk; None disables
        hetdimer_kmer_size    : k-mer size for heterodimer risk lookup
        gc_warning_buffer     : tolerance outside gc_range for flagging
        hairpin_dg_buffer     : ΔG flag buffer
        homodimer_dg_buffer   : ΔG flag buffer

    Returns:
        Returns total penalty and quality data dict for the given oligo
    """

    oligo_end   = pos + oligo_len
    overlap_seq = fragment[oligo_end - ov_len: oligo_end]
    min_tm, max_tm = tm_range

    # Reject a non positive length oligo or overlap
    if oligo_len <= 0 or ov_len <= 0 or len(overlap_seq) == 0:
        _warnings.warn(
            f"_evaluate_placement called with no oligo or no overlaps (pos={pos}, oligo_len={oligo_len}, ov_len={ov_len}, overlap_seq='{overlap_seq}') Returning max penalty to reject this placement",
            UserWarning,
            stacklevel=2,
        )
        return DP_PENALTY_FORCED, {
            'oligo_end':          oligo_end,
            'overlap_seq':        '',
            'tm':                 0.0,
            'oligo_seq':          '',
            'tm_ok':              False,
            'gc_ok':              False,
            'hp_ok':              True,
            'dn_ok':              True,
            'gc_warn':            ['Degenerate placement: zero-length overlap'],
            'hp_warn':            [],
            'dn_warn':            [],
            'vicinity_warn':      [],
            'hairpin_warn':       [],
            'homodimer_warn':     [],
            'gc_clamp_warn':      [],
            'body_hp_warn':       [],
            'body_dn_warn':       [],
            'misprime_hit_count': 0,
            'hetdimer_hit_count': 0,
            'penalty':            DP_PENALTY_FORCED,
        }

    # Use custom provided penalty values or fallback to default values if none provided
    p_tm          = penalties['tm']          if penalties else DP_PENALTY_TM
    p_gc          = penalties['gc']          if penalties else DP_PENALTY_GC
    p_homopolymer = penalties['homopolymer'] if penalties else DP_PENALTY_HOMOPOLYMER
    p_dinucleotide= penalties['dinucleotide']if penalties else DP_PENALTY_DINUCLEOTIDE
    p_hairpin     = penalties['hairpin']     if penalties else DP_PENALTY_HAIRPIN
    p_homodimer    = penalties['homodimer']   if penalties else DP_PENALTY_HOMODIMER
    p_gc_clamp     = penalties['gc_clamp']    if penalties else DP_PENALTY_GC_CLAMP
    p_misprime     = penalties['misprime']    if penalties else DP_PENALTY_MISPRIME
    p_hetdimer_map = (penalties.get('hetdimer_map', DP_PENALTY_HETDIMER_MAP)
                     if penalties else DP_PENALTY_HETDIMER_MAP)
    p_body_hp      = (penalties.get('body_homopolymer', DP_PENALTY_BODY_HOMOPOLYMER)
                     if penalties else DP_PENALTY_BODY_HOMOPOLYMER)
    p_body_dn      = (penalties.get('body_dinucleotide', DP_PENALTY_BODY_DINUCLEOTIDE)
                     if penalties else DP_PENALTY_BODY_DINUCLEOTIDE)
    sev_mult       = (penalties.get('severity_multiplier', DEFAULT_PENALTY_SEVERITY_MULTIPLIER)
                    if penalties else DEFAULT_PENALTY_SEVERITY_MULTIPLIER)

    # GC Check
    gc_val   = gc_percent(overlap_seq)
    gc_warn  = check_overlap_gc(overlap_seq, gc_range)
    gc_ok    = (gc_range[0] - gc_warning_buffer) <= gc_val <= (gc_range[1] + gc_warning_buffer)

    # Tm Check
    tm       = calc_tm(overlap_seq, mv_conc, dv_conc, dntp_conc, dna_conc)
    tm_ok    = (min_tm - tm_warning_buffer) <= tm <= (max_tm + tm_warning_buffer)

    # Homopolymer and dinucleotide Check
    hp_warn  = check_homopolymer(overlap_seq, max_homopolymer)
    dn_warn  = check_dinucleotide_repeat(overlap_seq, max_dinucleotide)
    hp_ok    = not hp_warn
    dn_ok    = not dn_warn

    penalty = 0.0

    if not tm_ok:
        effective_min = min_tm - tm_warning_buffer
        effective_max = max_tm + tm_warning_buffer
        tm_sev = (effective_min - tm) if tm < effective_min else (tm - effective_max)
        penalty += _scaled(p_tm, tm_sev, sev_mult)

    if gc_warn:  
        gc_sev = (gc_range[0] - gc_val) if gc_val < gc_range[0] else (gc_val - gc_range[1])
        penalty += _scaled(p_gc, max(0.0, gc_sev) / 10.0, sev_mult)  # normalise: 10% = severity 1

    if not hp_ok:
        penalty += _scaled(p_homopolymer, float(len(hp_warn)), sev_mult)

    if not dn_ok:
        penalty += _scaled(p_dinucleotide, float(len(dn_warn)), sev_mult)

    # Body homopolymer and dinucleotide check
    overlap_start_pos = oligo_end - ov_len
    oligo_body_seq    = fragment[pos : overlap_start_pos]
    body_hp_warn      = check_homopolymer(oligo_body_seq, max_homopolymer) if oligo_body_seq else []
    body_dn_warn      = check_dinucleotide_repeat(oligo_body_seq, max_dinucleotide) if oligo_body_seq else []

    if body_hp_warn:
        penalty += _scaled(p_body_hp, float(len(body_hp_warn)), sev_mult)
    if body_dn_warn:
        penalty += _scaled(p_body_dn, float(len(body_dn_warn)), sev_mult)

    # Homopolymer and dinucleotide vicinity check
    vicinity_warn = check_repeat_vicinity(
        fragment, overlap_start_pos, oligo_end,
        vicinity_buffer, max_homopolymer, max_dinucleotide,
        overlap_hp_results=hp_warn,
        overlap_dn_results=dn_warn,
    )

    # Misprime Check
    misprime_hit_count = 0
    if risk_map is not None:
        fwd_map, rev_map = risk_map
        kmer_size = misprime_kmer_size

        # Fragment-local safe zone (to avoid flagging intentional overlaps as a misprime)
        safe_lo = pos       - MISPRIME_SAFE_BUFFER
        safe_hi = oligo_end + MISPRIME_SAFE_BUFFER

        if direction == 'forward':
            # 3' end of ordered forward oligo = rightmost fragment-local pos
            map_idx = oligo_end - 1
            match_positions = fwd_map[map_idx] if 0 <= map_idx < len(fwd_map) else []
        else:
            # 3' end of ordered reverse oligo = leftmost top-strand pos = pos
            map_idx = pos
            match_positions = rev_map[map_idx] if 0 <= map_idx < len(rev_map) else []

        for site in match_positions:
            # off-target: k-mer span [site, site+kmer_size-1] does not overlap
            # the fragment-local safe zone [safe_lo, safe_hi)
            if site + kmer_size <= safe_lo or site >= safe_hi:
                misprime_hit_count += 1

        if misprime_hit_count > 0:
            # Severity = number of off-target hits (already a count)
            penalty += _scaled(p_misprime, float(misprime_hit_count - 1), sev_mult)

    # Step 5: Hairpin and Homodimer
    oligo_seq_top = fragment[pos: oligo_end]
    oligo_seq_ord = (oligo_seq_top if direction == 'forward'
                     else reverse_complement(oligo_seq_top))

    full_hp_warn, hp_sev = _hairpin_check_sliding(
        oligo_seq_ord, hairpin_dg_threshold,
        mv_conc, dv_conc, dntp_conc, dna_conc,
    )
    if hp_sev > 0.0:
        penalty += _scaled(p_hairpin, hp_sev, sev_mult)
    if hp_sev <= hairpin_dg_buffer:
        full_hp_warn = []

    full_hd_warn, hd_sev = _homodimer_check_sliding(
        oligo_seq_ord, homodimer_dg_threshold,
        mv_conc, dv_conc, dntp_conc, dna_conc,
    )
    if hd_sev > 0.0:
        penalty += _scaled(p_homodimer, hd_sev, sev_mult)
    if hd_sev <= homodimer_dg_buffer:
        full_hd_warn = []

    # GC clamp check
    gc_clamp_warn = check_gc_clamp(
        oligo_seq_ord, gc_clamp_min_gc, gc_clamp_length
    )
    if gc_clamp_warn:
        penalty += p_gc_clamp

    # Heterodimer check
    hetdimer_hit_count = 0
    if hetdimer_risk_map is not None:
        hd_safe_lo = pos       - HETDIMER_SAFE_BUFFER
        hd_safe_hi = oligo_end + HETDIMER_SAFE_BUFFER
        k          = hetdimer_kmer_size
        seen_sites: Set[int] = set()
        for j in range(len(oligo_seq_ord) - k + 1):
            rc_kmer = oligo_seq_ord[j:j + k].translate(_RC_TABLE)[::-1]
            for site in hetdimer_risk_map.get(rc_kmer, []):
                if site not in seen_sites and (
                    site + k <= hd_safe_lo or site >= hd_safe_hi
                ):
                    hetdimer_hit_count += 1
                    seen_sites.add(site)
        if hetdimer_hit_count > 0:
            penalty += _scaled(p_hetdimer_map, float(hetdimer_hit_count - 1), sev_mult)

    return penalty, {
        'oligo_end':            oligo_end,
        'overlap_seq':          overlap_seq,
        'tm':                   tm,
        'oligo_seq':            oligo_seq_ord,
        'tm_ok':                tm_ok,
        'gc_ok':                gc_ok,
        'hp_ok':                hp_ok,
        'dn_ok':                dn_ok,
        'gc_warn':              gc_warn,
        'hp_warn':              hp_warn,
        'dn_warn':              dn_warn,
        'vicinity_warn':        vicinity_warn,
        'hairpin_warn':         full_hp_warn,
        'homodimer_warn':       full_hd_warn,
        'gc_clamp_warn':        gc_clamp_warn,
        'body_hp_warn':         body_hp_warn,
        'body_dn_warn':         body_dn_warn,
        'misprime_hit_count':   misprime_hit_count,
        'hetdimer_hit_count':   hetdimer_hit_count,
        'penalty':              penalty,
    }

def _tile_fragment_dp(
    fragment: str, frag_idx: int,
    oligo_range: Tuple[int, int],
    overlap_range: Tuple[int, int],
    tm_range: Tuple[float, float],
    tm_warning_buffer: float,
    gc_range: Tuple[float, float],
    max_homopolymer: int,
    max_dinucleotide: int,
    vicinity_buffer: int,
    hairpin_dg_threshold: float,
    mv_conc: float, dv_conc: float, dntp_conc: float, dna_conc: float,
    frag_len_range: Optional[Tuple[int, int]],
    risk_map: Optional[Tuple[List[List[int]], List[List[int]]]] = None,
    frag_global_start: int = 0,
    penalties: Optional[dict] = None,
    homodimer_dg_threshold: float = DEFAULT_HOMODIMER_DG_THRESHOLD,
    gc_clamp_length: int = DEFAULT_GC_CLAMP_LENGTH,
    gc_clamp_min_gc: int = DEFAULT_GC_CLAMP_MIN_GC,
    misprime_kmer_size: int = DEFAULT_MISPRIME_KMER_SIZE,
    hetdimer_risk_map: Optional[Dict[str, List[int]]] = None,
    hetdimer_kmer_size: int = DEFAULT_HETDIMER_KMER_SIZE,
    gc_warning_buffer: float = DEFAULT_GC_WARNING_BUFFER,
    hairpin_dg_buffer: float = DEFAULT_HAIRPIN_DG_BUFFER,
    homodimer_dg_buffer: float = DEFAULT_HOMODIMER_DG_BUFFER,
    fwd_cassette: str = "",
    rev_cassette: str = "",
) -> List[dict]:
    """
    Tile one PCA fragment into alternating forward / reverse oligos using
    A* search to find the globally optimal placement with maximum speed.

    Graph model: the tiling is a DAG whose nodes are candidate oligo placements
    (start, length, strand) and whose edges are the valid overlap junctions between
    successive oligos. In the implementation, placements whose successor starts at
    the same fragment-local position on the same strand parity are merged into one
    state keyed by (pos, parity), so the search runs over those states from the single
    start state (0, 0) with non-negative oligo penalties as edge costs. Expanding a
    state evaluates every combination of oligo length and overlap length in the
    allowed ranges. Each state keeps the lowest accumulated penalty found so far
    (h = 0, Dijkstra-equivalent), and the cheapest terminal state is backtracked.

    Args:
        fragment          : DNA sequence of this fragment
        frag_idx          : 0-based fragment index (for record naming)
        oligo_range       : (min_bp, max_bp) for oligo length
        overlap_range     : (min_bp, max_bp) for overlap length
        tm_range          : (min_C, max_C) for overlap Tm
        gc_range          : (min_%, max_%) for overlap GC content
        max_homopolymer   : flag homopolymer runs >= this length
        max_dinucleotide  : flag dinucleotide repeats >= this count
        vicinity_buffer   : bp to scan around overlap boundary
        hairpin_dg_threshold : ΔG cutoff for hairpin flagging (kcal/mol)
        mv_conc, dv_conc, dntp_conc, dna_conc : reaction conditions
        frag_len_range    : (min_bp, max_bp) for fragment length warning
        risk_map          : pre-computed misprime risk scores from
                            _build_misprime_risk_map(). None = misprime
                            penalty disabled (backward-compatible default).
        frag_global_start : retained for API compatibility; unused (risk map is fragment-local)

    Returns:
        List of oligo record dicts (one per oligo) in 5'→3' order.
    """
    min_ol, max_ol = oligo_range
    min_ov, max_ov = overlap_range
    frag_len       = len(fragment)
    
    fwd_cassette_len = len(fwd_cassette)
    rev_cassette_len = len(rev_cassette)

    # Fragment-level size warning
    frag_size_warnings = []
    if frag_len_range:
        if frag_len < frag_len_range[0]:
            frag_size_warnings.append(
                f"Fragment length {frag_len} bp is below min ({frag_len_range[0]} bp)"
            )
        elif frag_len > frag_len_range[1]:
            frag_size_warnings.append(
                f"Fragment length {frag_len} bp is above max ({frag_len_range[1]} bp)"
            )
    # Iteration limit: scales with search space (200,000 floor for small fragments)
    _oligo_range_size   = max_ol - min_ol + 1
    _overlap_range_size = max_ov - min_ov + 1
    _max_iterations = max(200_000, frag_len * _oligo_range_size * _overlap_range_size * 5)

    # h = 0 everywhere (admissible: perfect sequence has zero penalty per oligo)
    # dp[(pos, parity)] = (g_cost, prev_pos, prev_parity, placement_data)
    INF = float('inf')
    dp: Dict[Tuple[int, int], Tuple[float, int, int, Optional[dict]]] = {}

    # Start state
    dp[(0, 0)] = (0.0, -1, -1, None)
    heap = [(0.0, 0.0, 0, 0)]   # (f=0.0, g=0.0, pos=0, parity=0)

    # Closed set: states whose optimal cost is finalised
    closed: set = set()

    iteration = 0

    while heap:
        iteration += 1

        # Runaway iteration guard
        if iteration > _max_iterations:
            raise RuntimeError(
                f"_tile_fragment_dp (A*) exceeded {_max_iterations:,} iterations on fragment {frag_idx+1} ({frag_len} bp), This indicates contradictory constraints or an unusually difficult sequence. Please review oligo_length_range, overlap_length_range, and overlap_tm_range, they may be mutually impossible for this sequence."
            )

        # Pop cheapest state
        f_val, g_cost, pos, parity = heapq.heappop(heap)

        state = (pos, parity)

        # Skip if already finalised
        if state in closed:
            continue
        closed.add(state)

        # Skip if stored cost is better than what we just popped
        stored_g = dp.get(state, (INF,))[0]
        if g_cost > stored_g + 1e-9:
            continue

        # Skip terminal bookkeeping entries re-popped from heap
        if pos >= frag_len:
            continue

        remaining = frag_len - pos
        direction = 'forward' if parity == 0 else 'reverse'
        
        max_ol_here = max_ol
        if pos == 0:
            max_ol_here -= fwd_cassette_len

        # Terminal state -> remaining fits in one oligo
        if remaining <= max_ol_here - rev_cassette_len:
            oligo_seq_top = fragment[pos:]

            # skip empty terminal slice (prevents primer3 C crash!!)
            if not oligo_seq_top:
                continue

            oligo_seq = (oligo_seq_top if direction == 'forward'
                         else reverse_complement(oligo_seq_top))

            final_oligo_len = len(oligo_seq_top) + rev_cassette_len
            if pos == 0:
                final_oligo_len += fwd_cassette_len
                
            term_warnings = list(frag_size_warnings) if pos == 0 else []
            if final_oligo_len < min_ol:
                term_warnings.append(
                    f"Last oligo {final_oligo_len} bp < min {min_ol} bp"
                )
            hp_w, _    = _hairpin_check_sliding(
                oligo_seq, hairpin_dg_threshold,
                mv_conc, dv_conc, dntp_conc, dna_conc,
            )
            hd_w, _    = _homodimer_check_sliding(
                oligo_seq, homodimer_dg_threshold,
                mv_conc, dv_conc, dntp_conc, dna_conc,
            )
            gc_clamp_w = check_gc_clamp(oligo_seq, gc_clamp_min_gc, gc_clamp_length)

            # For the terminal oligo there is no overlap, so the entire oligo is
            # the body. We compute body flags here so they are not silently False.
            term_body_hp_warn = check_homopolymer(oligo_seq, max_homopolymer)
            term_body_dn_warn = check_dinucleotide_repeat(oligo_seq, max_dinucleotide)

            terminal_data = {
                'oligo_seq':      oligo_seq,
                'frag_end':       frag_len,
                'warnings':       term_warnings + hp_w + hd_w + gc_clamp_w,
                'is_last':        True,
                'hairpin_warn':   hp_w,
                'homodimer_warn': hd_w,
                'gc_clamp_warn':  gc_clamp_w,
                'tm_ok': True, 'gc_ok': True, 'hp_ok': True, 'dn_ok': True,
                'overlap_seq': None, 'tm': None,
                'gc_warn': [], 'hp_warn': [], 'dn_warn': [], 'vicinity_warn': [],
                'body_hp_warn':   term_body_hp_warn,
                'body_dn_warn':   term_body_dn_warn,
            }
            terminal_state = (frag_len, 1 - parity)
            old_g = dp.get(terminal_state, (INF,))[0]
            if g_cost < old_g:
                dp[terminal_state] = (g_cost, pos, parity, terminal_data)
                heapq.heappush(heap, (g_cost, g_cost, frag_len, 1 - parity))
            continue

        # Enumerate all valid placements from this state
        for oligo_len in range(min_ol, min(max_ol_here, remaining - 1) + 1):
            oligo_end   = pos + oligo_len
            max_ov_here = min(max_ov, oligo_len - 1)

            for ov_len in range(min_ov, max_ov_here + 1):
                penalty, ev = _evaluate_placement(
                    fragment, pos, oligo_len, ov_len, direction,
                    tm_range, tm_warning_buffer, gc_range, max_homopolymer, max_dinucleotide,
                    vicinity_buffer, hairpin_dg_threshold,
                    mv_conc, dv_conc, dntp_conc, dna_conc,
                    risk_map=risk_map,
                    frag_global_start=frag_global_start,
                    penalties=penalties,
                    homodimer_dg_threshold=homodimer_dg_threshold,
                    gc_clamp_length=gc_clamp_length,
                    gc_clamp_min_gc=gc_clamp_min_gc,
                    misprime_kmer_size=misprime_kmer_size,
                    hetdimer_risk_map=hetdimer_risk_map,
                    hetdimer_kmer_size=hetdimer_kmer_size,
                    gc_warning_buffer=gc_warning_buffer,
                    hairpin_dg_buffer=hairpin_dg_buffer,
                    homodimer_dg_buffer=homodimer_dg_buffer,
                )

                next_pos    = oligo_end - ov_len
                next_pos    = max(next_pos, pos + 1)
                next_parity = 1 - parity
                new_g       = g_cost + penalty

                # Unreachable given max() above; guards against future refactors
                if next_pos <= pos:
                    continue

                next_state = (next_pos, next_parity)
                old_g      = dp.get(next_state, (INF,))[0]

                if new_g < old_g:
                    # Build warnings for this placement
                    w = list(frag_size_warnings) if pos == 0 else []
                    w += ev['vicinity_warn'] + ev['hairpin_warn'] + ev['homodimer_warn'] + ev['gc_clamp_warn']
                    if not ev['tm_ok']:
                        bound = tm_range[0] if ev['tm'] < (tm_range[0] - tm_warning_buffer) else tm_range[1]
                        direction_word = ('BELOW minimum' if ev['tm'] < (tm_range[0] - tm_warning_buffer)
                                          else 'ABOVE maximum')
                        w.append(f"Overlap Tm {ev['tm']:.1f}C is "
                                 f"{direction_word} {bound}C")
                    w += ev['gc_warn'] + ev['hp_warn'] + ev['dn_warn']

                    placement_data = {
                        'oligo_seq':      ev['oligo_seq'],
                        'frag_end':       oligo_end,
                        'overlap_seq':    ev['overlap_seq'],
                        'tm':             ev['tm'],
                        'warnings':       w,
                        'is_last':        False,
                        'hairpin_warn':   ev['hairpin_warn'],
                        'homodimer_warn': ev['homodimer_warn'],
                        'gc_clamp_warn':  ev['gc_clamp_warn'],
                        'tm_ok':          ev['tm_ok'],
                        'gc_ok':          ev['gc_ok'],
                        'hp_ok':          ev['hp_ok'],
                        'dn_ok':          ev['dn_ok'],
                        'gc_warn':        ev['gc_warn'],
                        'hp_warn':        ev['hp_warn'],
                        'dn_warn':        ev['dn_warn'],
                        'vicinity_warn':  ev['vicinity_warn'],
                        'body_hp_warn':   ev['body_hp_warn'],
                        'body_dn_warn':   ev['body_dn_warn'],
                    }
                    dp[next_state] = (new_g, pos, parity, placement_data)

                    # Push to heap with f = g + h; since h = 0, f = g
                    heapq.heappush(heap, (new_g, new_g, next_pos, next_parity))

    # Verify at least one terminal state was reached
    terminal_states = [(frag_len, p) for p in (0, 1)
                       if dp.get((frag_len, p), (INF,))[0] < INF]

    if not terminal_states:
        raise RuntimeError(
            f"_tile_fragment_dp (A*) could not tile fragment {frag_idx+1} "
            f"(length {frag_len} bp) within the given constraints. "
            f"No valid tiling path was found. "
            f"Please review: oligo_length_range={oligo_range}, "
            f"overlap_length_range={overlap_range}, "
            f"overlap_tm_range={tm_range}. "
            f"These constraints may be mutually impossible for this sequence."
        )

    # Select the terminal state with lowest g cost
    best_terminal = min(terminal_states, key=lambda s: dp[s][0])

    # Backtrack to reconstruct the optimal oligo sequence
    path = []
    state = best_terminal
    while state != (-1, -1):
        cost, prev_pos, prev_parity, data = dp[state]
        if data is not None:
            path.append((state, data))
        state = (prev_pos, prev_parity)
    path.reverse()

    # Build records from the reconstructed path
    oligos    = []
    oligo_idx = 0

    for (pos, parity), data in path:
        direction = 'forward' if oligo_idx % 2 == 0 else 'reverse'

        if oligo_idx == 0:
            frag_start_here = 0
            prev_g = 0.0
        else:
            prev_state, _ = path[oligo_idx - 1]
            frag_start_here = prev_state[0]
            prev_g = dp[prev_state[0], prev_state[1]][0]

        frag_end_here = data['frag_end']
        state         = (pos, parity)
        oligo_penalty = round(dp[state][0] - prev_g, 6)

        rec = _make_record(
            frag_idx, oligo_idx, direction,
            frag_start_here, frag_end_here,
            data['oligo_seq'],
            data['overlap_seq'],
            round(data['tm'], 2) if data['tm'] is not None else None,
            data['warnings'],
            is_last=data['is_last'],
            tm_flag=not data['tm_ok'],
            gc_flag=not data['gc_ok'],
            homopolymer_flag=not data['hp_ok'],
            dinucleotide_flag=not data['dn_ok'],
            hairpin_flag=bool(data['hairpin_warn']),
            homodimer_flag=bool(data['homodimer_warn']),
            gc_clamp_flag=bool(data['gc_clamp_warn']),
            body_homopolymer_flag=bool(data.get('body_hp_warn')),
            body_dinucleotide_flag=bool(data.get('body_dn_warn')),
        )
        rec['penalty_score'] = oligo_penalty
        oligos.append(rec)
        oligo_idx += 1

    return oligos

# SECTION 6b: PARALLEL TILING WORKER
# Falls back to sequential execution if multiprocessing is unavailable

def _tile_fragment_worker(args: tuple) -> Tuple[int, List[dict]]:
    """
    ProcessPoolExecutor worker: unpacks args tuple and delegates to _tile_fragment_dp.
    Must be at module level (picklable). Returns (frag_idx, records).
    """
    (fragment, frag_idx, oligo_range, overlap_range, tm_range, tm_warning_buffer, gc_range,
     max_homopolymer, max_dinucleotide, vicinity_buffer, hairpin_dg_threshold,
     mv_conc, dv_conc, dntp_conc, dna_conc, frag_len_range,
     risk_map, frag_global_start,
     penalties, homodimer_dg_threshold, gc_clamp_length, gc_clamp_min_gc,
     misprime_kmer_size, hetdimer_risk_map, hetdimer_kmer_size,
     gc_warning_buffer, hairpin_dg_buffer, homodimer_dg_buffer,
     fwd_cassette, rev_cassette) = args

    records = _tile_fragment_dp(
        fragment=fragment, frag_idx=frag_idx,
        oligo_range=oligo_range, overlap_range=overlap_range,
        tm_range=tm_range, tm_warning_buffer=tm_warning_buffer, gc_range=gc_range,
        max_homopolymer=max_homopolymer, max_dinucleotide=max_dinucleotide,
        vicinity_buffer=vicinity_buffer, hairpin_dg_threshold=hairpin_dg_threshold,
        mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=dna_conc,
        frag_len_range=frag_len_range,
        risk_map=risk_map,
        frag_global_start=frag_global_start,
        penalties=penalties,
        homodimer_dg_threshold=homodimer_dg_threshold,
        gc_clamp_length=gc_clamp_length,
        gc_clamp_min_gc=gc_clamp_min_gc,
        misprime_kmer_size=misprime_kmer_size,
        hetdimer_risk_map=hetdimer_risk_map,
        hetdimer_kmer_size=hetdimer_kmer_size,
        gc_warning_buffer=gc_warning_buffer,
        hairpin_dg_buffer=hairpin_dg_buffer,
        homodimer_dg_buffer=homodimer_dg_buffer,
        fwd_cassette=fwd_cassette,
        rev_cassette=rev_cassette,
    )
    return frag_idx, records

# SECTION 7a: LOGGING UTILITIES

def _fmt_elapsed(seconds: float) -> str:
    """Formats a seconds duration to a human-readable string"""
    m = int(seconds // 60)
    s = int(seconds % 60)
    return f"{m}m {s:02d}s" if m > 0 else f"{s}s"

# SECTION 7b: ASSEMBLY VALIDATION

def _validate_fragment_assembly(fragment_seq: str, fragment_records: list) -> None:
    """
    Simulates PCA for a single fragment and confirms the designed oligos perfectly
    reconstruct the intended fragment sequence. Raises ValueError on failure.
    """
    if not fragment_records:
        return
        
    frag_num = fragment_records[0]['fragment']
    
    # Sort oligos by their position in the fragment
    sorted_records = sorted(fragment_records, key=lambda r: r['oligo_in_frag'])
    
    reconstructed_seq = ""
    for i, rec in enumerate(sorted_records):
        oligo_seq = str(rec['oligo_sequence']).upper()
        direction = rec['direction']
        
        if direction == 'forward':
            top_strand_seq = oligo_seq
        else:
            top_strand_seq = reverse_complement(oligo_seq).upper()
            
        if i == 0:
            reconstructed_seq += top_strand_seq
        else:
            prev_overlap = str(sorted_records[i-1]['overlap_sequence']).upper()
            if not reconstructed_seq.endswith(prev_overlap):
                raise ValueError(f"Fragment {frag_num} validation failed: Overlap mismatch for oligo {rec['oligo_name']}.")
            if not top_strand_seq.startswith(prev_overlap):
                raise ValueError(f"Fragment {frag_num} validation failed: Oligo {rec['oligo_name']} does not start with expected overlap.")
            
            reconstructed_seq += top_strand_seq[len(prev_overlap):]

    if reconstructed_seq != fragment_seq.upper():
         raise ValueError(f"Fragment {frag_num} validation failed: Reconstructed sequence does not match intended fragment sequence.")

def _validate_full_assembly(
    original_seq: str,
    fragment_df: pd.DataFrame,
    assembly_method: str,
    junction_overlap_length: int,
    gg_overhangs: Optional[List[str]] = None,
    goldengate_enzyme: str = "BsaI"
) -> None:
    """
    Simulates the final assembly reaction of all produced fragments back into the original full-length sequence.
    Raises ValueError on failure.
    """
    if fragment_df.empty:
        return
        
    fragments = fragment_df['fragment_sequence'].tolist()
    
    if assembly_method in ["stitching_pcr", "gibson"]:
        reconstructed_seq = fragments[0].upper()
        for i in range(1, len(fragments)):
            curr_frag = fragments[i].upper()
            overlap = curr_frag[:junction_overlap_length]
            if not reconstructed_seq.endswith(overlap):
                raise ValueError(f"Full Assembly Validation failed: Junction overlap mismatch between fragment {i} and {i+1}.")
            reconstructed_seq += curr_frag[junction_overlap_length:]
            
        if reconstructed_seq != original_seq.upper():
            raise ValueError("Full Assembly Validation failed: Reassembled sequence does not match original sequence.")
            
    elif assembly_method == "golden_gate":
        enzyme = GOLDEN_GATE_ENZYMES[goldengate_enzyme]
        rec_seq = enzyme["recognition"].upper()
        rc_rec = enzyme["rc_recognition"].upper()
        spacer = enzyme["fwd_spacer"].upper()
        extras = GOLDEN_GATE_EXTRA_BASES.upper()
        
        digested_frags = []
        for i, frag in enumerate(fragments):
            frag_up = frag.upper()
            left_ov = gg_overhangs[i]
            right_ov = gg_overhangs[i+1]
            
            left_strip = 0
            if left_ov is not None:
                left_ov = left_ov.upper()
                left_cassette = extras + rec_seq + spacer
                if not frag_up.startswith(left_cassette):
                    raise ValueError(f"Full Assembly Validation failed: Fragment {i+1} missing or malformed 5' Golden Gate cassette.")
                left_strip = len(left_cassette)
                
            right_strip = 0
            if right_ov is not None:
                right_ov = right_ov.upper()
                right_cassette = spacer + rc_rec + extras
                if not frag_up.endswith(right_cassette):
                    raise ValueError(f"Full Assembly Validation failed: Fragment {i+1} missing or malformed 3' Golden Gate cassette.")
                right_strip = len(right_cassette)
                
            digested_top = frag_up[left_strip : len(frag_up) - right_strip]
            digested_frags.append(digested_top)
            
        reconstructed_seq = digested_frags[0]
        for i in range(1, len(digested_frags)):
            prev_frag = reconstructed_seq
            curr_frag = digested_frags[i]
            overlap = gg_overhangs[i].upper()
            
            if not prev_frag.endswith(overlap):
                raise ValueError(f"Full Assembly Validation failed: Fragment {i} does not end with expected 4-bp overhang {overlap}.")
            if not curr_frag.startswith(overlap):
                raise ValueError(f"Full Assembly Validation failed: Fragment {i+1} does not start with expected 4-bp overhang {overlap}.")
                
            reconstructed_seq += curr_frag[len(overlap):]
            
        # Strip terminal overhangs for comparison with original sequence
        test_seq = reconstructed_seq
        if gg_overhangs[0] is not None:
            first_ov = gg_overhangs[0].upper()
            if test_seq.startswith(first_ov):
                test_seq = test_seq[len(first_ov):]
        if gg_overhangs[-1] is not None:
            last_ov = gg_overhangs[-1].upper()
            if test_seq.endswith(last_ov):
                test_seq = test_seq[:-len(last_ov)]
            
        if test_seq != original_seq.upper():
            raise ValueError(f"Full Assembly Validation failed: Reassembled Golden Gate sequence does not match original sequence. Length reconstructed: {len(test_seq)}, original: {len(original_seq)}.")


# SECTION 7c: MAIN PUBLIC FUNCTION

def design_oligos(
    sequence: str,
    # Preset mode
    mode:                         str                       = "custom",
    # Structural parameters
    fragment_length_range:    Optional[Tuple[int, int]]  = DEFAULT_FRAGMENT_LENGTH_RANGE,
    oligo_length_range:       Tuple[int, int]            = DEFAULT_OLIGO_LENGTH_RANGE,
    overlap_length_range:     Tuple[int, int]            = DEFAULT_OVERLAP_LENGTH_RANGE,
    # Quality thresholds
    overlap_tm_range:         Tuple[float, float]        = DEFAULT_OVERLAP_TM_RANGE,
    tm_warning_buffer:        float                      = DEFAULT_TM_WARNING_BUFFER,
    gc_warning_buffer:        float                      = DEFAULT_GC_WARNING_BUFFER,
    overlap_gc_range:         Tuple[float, float]        = DEFAULT_OVERLAP_GC_RANGE,
    max_homopolymer_run:      int                        = DEFAULT_MAX_HOMOPOLYMER_RUN,
    max_dinucleotide_repeat:  int                        = DEFAULT_MAX_DINUCLEOTIDE_REPEAT,
    repeat_vicinity_buffer:   int                        = DEFAULT_REPEAT_VICINITY_BUFFER,
    hairpin_dg_threshold:     float                      = DEFAULT_HAIRPIN_DG_THRESHOLD,
    hairpin_dg_buffer:        float                      = DEFAULT_HAIRPIN_DG_BUFFER,
    homodimer_dg_threshold:   float                      = DEFAULT_HOMODIMER_DG_THRESHOLD,
    homodimer_dg_buffer:      float                      = DEFAULT_HOMODIMER_DG_BUFFER,
    gc_clamp_length:          int                        = DEFAULT_GC_CLAMP_LENGTH,
    gc_clamp_min_gc:          int                        = DEFAULT_GC_CLAMP_MIN_GC,
    heterodimer_dg_threshold: float                      = DEFAULT_HETERODIMER_DG_THRESHOLD,
    heterodimer_dg_buffer:    float                      = DEFAULT_HETERODIMER_DG_BUFFER,
    heterodimer_neighbor_range: int                      = DEFAULT_HETERODIMER_NEIGHBOR_RANGE,
    # Misprime detection
    misprime_length:          int                        = DEFAULT_MISPRIME_LENGTH,
    misprime_tip_length:      int                        = DEFAULT_MISPRIME_TIP_LENGTH,
    misprime_max_mismatches:  int                        = DEFAULT_MISPRIME_MAX_MISMATCHES,
    misprime_kmer_size:       int                        = DEFAULT_MISPRIME_KMER_SIZE,
    # DP penalty scores
    penalty_tm:               float                      = DP_PENALTY_TM,
    penalty_gc:               float                      = DP_PENALTY_GC,
    penalty_homopolymer:      float                      = DP_PENALTY_HOMOPOLYMER,
    penalty_dinucleotide:     float                      = DP_PENALTY_DINUCLEOTIDE,
    penalty_hairpin:          float                      = DP_PENALTY_HAIRPIN,
    penalty_homodimer:        float                      = DP_PENALTY_HOMODIMER,
    penalty_gc_clamp:         float                      = DP_PENALTY_GC_CLAMP,
    penalty_misprime:         float                      = DP_PENALTY_MISPRIME,
    penalty_hetdimer_map:     float                      = DP_PENALTY_HETDIMER_MAP,
    penalty_body_homopolymer: float                      = DP_PENALTY_BODY_HOMOPOLYMER,
    penalty_body_dinucleotide:float                      = DP_PENALTY_BODY_DINUCLEOTIDE,
    penalty_severity_multiplier: float                   = DEFAULT_PENALTY_SEVERITY_MULTIPLIER,
    # Heterodimer RC k-mer map
    hetdimer_kmer_size:       int                        = DEFAULT_HETDIMER_KMER_SIZE,
    # Fragment cut scoring
    cut_score_window:         int                        = 20,
    cut_score_repeat_penalty: float                      = 1000.0,
    # Thermodynamic reaction conditions
    mv_conc:                  float                      = DEFAULT_MV_CONC,
    dv_conc:                  float                      = DEFAULT_DV_CONC,
    dntp_conc:                float                      = DEFAULT_DNTP_CONC,
    dna_conc:                 float                      = DEFAULT_DNA_CONC,
    # Output
    output_csv:               Optional[str]              = DEFAULT_OUTPUT_CSV,
    output_log:               Optional[str]              = None,
    # Inter-fragment junction
    junction_overlap_length:  int                        = DEFAULT_JUNCTION_OVERLAP_LENGTH,
    junction_tm_range:        Tuple[float, float]        = DEFAULT_JUNCTION_TM_RANGE,
    # Fragment boundary search
    cut_search_window:        int                        = DEFAULT_CUT_SEARCH_WINDOW,
    # Outer PCR primers
    outer_primer_length_range: Tuple[int, int]           = DEFAULT_OUTER_PRIMER_LENGTH_RANGE,
    outer_primer_tm_range:     Tuple[float, float]       = DEFAULT_OUTER_PRIMER_TM_RANGE,
    # Assembly method
    assembly_method:           str                       = DEFAULT_ASSEMBLY_METHOD,
    # Golden Gate Assembly
    goldengate_enzyme:          str                       = DEFAULT_GOLDENGATE_ENZYME,
    # Workflow control
    run_oligos:               bool                       = True,
    run_primers:              bool                       = True,
    fragment_indices:         Optional[List[int]]        = None,
    # Parallelisation
    n_workers:                Optional[int]              = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Design PCA oligos, junction overlaps, and outer PCR primers for a DNA sequence.

    Parameters
    ----------
    sequence : str
        Input DNA sequence (5'→3', top strand). Accepts raw ACGT strings, single-entry FASTA,
        and multi-line GenBank format. Whitespace and digits stripped automatically.
        Only A, T, G, C supported. Multi-entry FASTA: only first entry used (UserWarning emitted).

    mode : str
        Preset mode that controls structural parameter behaviour. One of:
          "custom"       : (default) all parameters are user-controlled.
          "short_oligos" : oligo 40–60 bp, overlap 15–25 bp, fragment 400–500 bp.
          "long_oligos"  : oligo 680–700 bp, overlap 30–40 bp, fragment 3000–4000 bp.
        Presets override only the three structural ranges (oligo length,
        overlap length, fragment length). All other params remain
        user-controllable or at their function-signature defaults.

    fragment_length_range : tuple (min_bp, max_bp) or None
        Target fragment length. None = treat entire sequence as one fragment.
        Cut points are sequence-aware (avoids homopolymers / GC extremes). Default: (400, 500)

    oligo_length_range : tuple (min_bp, max_bp)
        Oligo length range. Shorter oligos tried first (lower synthesis cost). Default: (40, 60)

    overlap_length_range : tuple (min_bp, max_bp)
        Overlap region shared between consecutive oligos (hybridises during PCA). Default: (15, 25)

    overlap_tm_range : tuple (min_°C, max_°C)
        Tm range for overlap regions (SantaLucia 1998 nearest-neighbor, salt-corrected). Default: (60.0, 65.0)

    tm_warning_buffer : float
        Tolerance in °C added on each side of overlap_tm_range before an overlap Tm is flagged. Default: 1.0

    overlap_gc_range : tuple (min_%, max_%)
        GC content range for overlap regions. Default: (30.0, 65.0)

    gc_warning_buffer : float
        Tolerance in percentage points added on each side of overlap_gc_range before an overlap GC% is flagged. Default: 2.0

    max_homopolymer_run : int
        Flag homopolymer runs >= this length in overlap (penalised in A*) and oligo body
        (softer penalty, reported as flag_body_homopolymer). Default: 4

    max_dinucleotide_repeat : int
        Flag dinucleotide repeats >= this many repetitions in the overlap. Default: 3

    repeat_vicinity_buffer : int
        Extra bp scanned on each side of the overlap for repeats just outside it. Produces a proximity
        warning only (no flag, no penalty). Default: 5

    hairpin_dg_threshold : float
        Flag oligos with intramolecular hairpin dG <= this value (kcal/mol). Default: -9.0

    hairpin_dg_buffer : float
        Extra kcal/mol beyond hairpin_dg_threshold that is tolerated before a hairpin is flagged. Default: 1.0

    homodimer_dg_buffer : float
        Extra kcal/mol beyond homodimer_dg_threshold that is tolerated before a homodimer is flagged. Default: 1.0

    heterodimer_dg_buffer : float
        Extra kcal/mol beyond heterodimer_dg_threshold that is tolerated before a heterodimer is flagged. Default: 1.0

    misprime_length : int
        Length of the 3' probe used for post-tiling misprime scan. Default: 20

    misprime_tip_length : int
        Exact-match tip length at probe 3' end. Sites must match this tip before
        the rest of the probe is evaluated. Default: 8

    misprime_max_mismatches : int
        Max mismatches in the non-tip probe region for a site to count as a misprime. Default: 2

    mv_conc : float
        Monovalent salt concentration in mM (Na+/K+) for Tm calculations. Default: 50.0

    dv_conc : float
        Divalent salt concentration in mM (Mg2+). Default: 3.0

    dntp_conc : float
        dNTP concentration in mM (chelates Mg2+, reducing effective free Mg2+). Default: 0.8

    dna_conc : float
        Total oligo strand concentration in nM. Default: 200.0

    output_csv : str or None
        Output file path. Also writes <path>_fragments.csv and <path>_outer_primers.csv.
        None = skip saving. Default: "mySequence.csv"

    junction_overlap_length : int
        Shared overlap between consecutive fragments for stitching PCR
        (last N bp of frag N == first N bp of frag N+1). Default: 25

    junction_tm_range : tuple (min_°C, max_°C)
        Tm range for junction overlaps (wider than overlap_tm_range). Default: (55.0, 68.0)

    homodimer_dg_threshold : float
        Flag oligos with homodimer dG <= this value (kcal/mol). Default: -9.0

    gc_clamp_length : int
        Number of 3' terminal bases examined for GC clamp check. Default: 3

    gc_clamp_min_gc : int
        Minimum G/C count required in the last gc_clamp_length bases. Default: 1

    heterodimer_dg_threshold : float
        Flag oligo pairs with heterodimer dG <= this value (kcal/mol).
        ±1 neighbors skipped (intentional overlap partners). Default: -9.0

    heterodimer_neighbor_range : int
        Neighbor range for heterodimer checks on each side (±1 always skipped). Default: 2

    misprime_kmer_size : int
        k-mer seed length for the pre-tiling per-fragment misprime risk map. Default: 12

    penalty_tm : float
        A* penalty when overlap Tm is outside overlap_tm_range. Default: 20

    penalty_gc : float
        A* penalty when overlap GC% is outside overlap_gc_range. Default: 8

    penalty_homopolymer : float
        A* penalty for a homopolymer run in the overlap. Default: 10

    penalty_dinucleotide : float
        A* penalty for a dinucleotide repeat in the overlap. Default: 6

    penalty_hairpin : float
        A* penalty when the full oligo has a stable hairpin. Default: 10

    penalty_homodimer : float
        A* penalty when the full oligo has homodimer risk. Default: 6

    penalty_gc_clamp : float
        A* penalty for weak 3' GC content. Default: 5

    penalty_misprime : float
        A* penalty per off-target k-mer hit in the pre-tiling misprime risk map. Default: 6

    hetdimer_kmer_size : int
        k-mer seed length for the pre-tiling heterodimer RC risk map. Default: 12

    penalty_hetdimer_map : float
        A* penalty per unique RC k-mer hit in the pre-tiling hetdimer risk map. Default: 8

    penalty_body_homopolymer : float
        Soft A* steering penalty for homopolymer runs in the oligo body (non-overlap).
        Flagged as flag_body_homopolymer. Set to 0 to disable steering. Default: 3

    penalty_body_dinucleotide : float
        Same as penalty_body_homopolymer but for dinucleotide repeats in the body.
        Flagged as flag_body_dinucleotide. Default: 3

    penalty_severity_multiplier : float
        Scales penalties proportionally to how far a threshold is violated.
        Formula: penalty = base * (1 + severity_multiplier * severity)
        0.0 = flat; 1.0 = proportional (default); 2.0 = steep. Default: 1.0

    cut_search_window : int
        bp searched on either side of the target cut position for the lowest-penalty fragment boundary. Default: 30

    cut_score_window : int
        Sequence window (bp) examined around each candidate fragment cut point. Default: 20

    cut_score_repeat_penalty : float
        Penalty per repeat element in the cut_score_window. Default: 1000.0

    run_oligos : bool
        False = skip tiling, misprime, heterodimer, junction, and primer design.
        Returns a fragment-only table. Default: True

    run_primers : bool
        False = skip outer PCR primer design. Oligo and junction tables still produced. Default: True

    fragment_indices : list of int or None
        1-based list of fragments to tile. None = tile all. Default: None

    outer_primer_length_range : tuple (min_bp, max_bp)
        Outer PCR primer length range. For each fragment end, the length whose Tm is closest to the middle of
        outer_primer_tm_range is chosen. Default: (15, 35)

    outer_primer_tm_range : tuple (min_°C, max_°C)
        Tm range for outer PCR primers. Default: (55.0, 68.0)

    assembly_method : str
        Method for joining PCA fragments after outer PCR. Default: "gibson"

        "stitching_pcr"
            Fragments share junction_overlap_length bp for overlap-extension PCR.

        "gibson"
            Gibson Assembly (50°C isothermal). Auto-enforces: (1) junction_overlap_length
            raised to 40 bp minimum; (2) junction_tm_range minimum relaxed to 48°C.

        "golden_gate"
            Golden Gate (Type IIS enzymes, 4-bp sticky ends). Enzyme auto-switched if the
            preferred enzyme has internal recognition sites. ValueError if all three blocked.

    goldengate_enzyme : str
        Type IIS enzyme for Golden Gate: "BsaI" | "BsmBI" | "BbsI". Auto-switched if the
        chosen enzyme has internal recognition sites in the input. Default: "BsaI"

    output_log : str or None
        Log file path. Captures everything printed by the run (not Python warnings). Always written before returning,
        including on error ([ERROR] line appended). Auto-derived from output_csv if None.
        Set both to None to suppress disk output. Default: None

    n_workers : int or None
        Parallel worker processes for fragment tiling. None = all cores. 1 = sequential. Default: None

    Returns
    -------
    tuple : (oligo_df, fragment_df, outer_primer_df)

    oligo_df : pd.DataFrame
        One row per oligo. Columns: global_oligo_num, oligo_name, fragment, oligo_in_frag,
        direction, frag_start, frag_end, global_start, global_end, oligo_length, oligo_gc_%,
        oligo_sequence, overlap_sequence, overlap_length, overlap_tm_C, overlap_gc_%,
        flag_tm_out_of_range, flag_overlap_gc, flag_homopolymer, flag_dinucleotide_repeat,
        flag_hairpin, flag_homodimer, flag_gc_clamp, flag_heterodimer, flag_misprime,
        flag_body_homopolymer, flag_body_dinucleotide, penalty_score, is_last_in_frag, warnings.
        penalty_score reflects the A* DP edge cost at tiling time (Tm, GC, hairpin, homodimer,
        misprime risk map, heterodimer risk map, GC clamp, body repeats). Post-tiling checks
        (flag_misprime probe scan, flag_heterodimer dG check) do not contribute to this score.

    fragment_df : pd.DataFrame
        One row per fragment. Always returned. Contains fragment boundaries, junction columns
        (assembly-method-specific), and per-fragment statistics.

    outer_primer_df : pd.DataFrame
        One row per primer pair per fragment. Empty DataFrame if run_primers=False.

    .. note::
        output_log is always written to disk before returning or raising. Read it via:
            log_text = open(output_log, encoding='utf-8').read()
    """
    # Start the timer
    _t_start  = time.time()

    # ── Preset mode enforcement ──────────────────────────────────────
    if mode not in VALID_MODES:
        raise ValueError(
            f"Invalid mode '{mode}'. Must be one of: {', '.join(VALID_MODES)}"
        )

    if mode != "custom":
        preset = _PRESETS[mode]
        # Apply structural overrides from preset
        fragment_length_range       = preset["fragment_length_range"]
        oligo_length_range          = preset["oligo_length_range"]
        overlap_length_range        = preset["overlap_length_range"]
        # All other params (quality, penalties, thermo, assembly, etc.)
        # remain as passed by the caller or their function-signature defaults.

    # Log capture: everything printed through _print() below goes to the console and to an
    # in-memory buffer that is saved as the run log (Python warnings are shown on the console only).
    import io as _io, re as _re
    if output_log is None and output_csv is not None:
        _log_stem  = _re.sub(r'\.csv$', '', output_csv, flags=_re.IGNORECASE)
        output_log = _log_stem + '_log.txt'
    _log_buf     = _io.StringIO()
    _real_stdout = sys.stdout  # captured once at call time; never overwritten

    def _print(*args, **kwargs) -> None:
        """Thread-safe drop-in for print(): writes to log buffer AND real stdout."""
        import io as _io2
        _tmp = _io2.StringIO()
        kwargs_copy = {k: v for k, v in kwargs.items() if k != 'file'}
        print(*args, file=_tmp, **kwargs_copy)
        line = _tmp.getvalue()
        _log_buf.write(line)
        _real_stdout.write(line)
        _real_stdout.flush()

    _log_written = [False]   # guard against double-write

    def _write_log() -> None:
        """Flush log buffer to disk (thread-safe: no global stdout mutation)"""
        if _log_written[0]:
            return
        _log_written[0] = True
        if output_log:
            try:
                with open(output_log, 'w', encoding='utf-8') as _lf:
                    _lf.write(_log_buf.getvalue())
                print(f'  Run log saved to      : {output_log}')
            except OSError as _e:
                print(f'  [WARNING] Could not write log file {output_log}: {_e}')

    class _LogGuard:
        """guarantees _write_log() is called on any exit mode (return, raise, or crash)"""
        def __del__(self):
            _write_log()
    _log_guard = _LogGuard()  # lives until design_oligos() returns or raises

    def _raise(exc: Exception) -> None:
        """Write [ERROR] to log, flush to disk, then raise"""
        _log_buf.write(f'\n[ERROR] {type(exc).__name__}: {exc}\n')
        _write_log()
        raise exc

    # Full parameter log
    _print("=" * 60)
    _print("  OLIGOFIT: Full Parameter Log")
    _print("=" * 60)
    _print(f"  mode                      : {mode}")
    _print("  -- Structural -----------------------------------------------")
    _print(f"  sequence_length           : {len(sequence)} bp")
    _print(f"  fragment_length_range     : {fragment_length_range}")
    _print(f"  oligo_length_range        : {oligo_length_range}")
    _print(f"  overlap_length_range      : {overlap_length_range}")
    _print("  -- Quality thresholds ----------------------------------------")
    _print(f"  overlap_tm_range          : {overlap_tm_range}")
    _print(f"  tm_warning_buffer         : {tm_warning_buffer} C")
    _print(f"  gc_warning_buffer         : {gc_warning_buffer} %")
    _print(f"  overlap_gc_range          : {overlap_gc_range}")
    _print(f"  max_homopolymer_run       : {max_homopolymer_run}")
    _print(f"  max_dinucleotide_repeat   : {max_dinucleotide_repeat}")
    _print(f"  repeat_vicinity_buffer    : {repeat_vicinity_buffer}")
    _print(f"  hairpin_dg_threshold      : {hairpin_dg_threshold} kcal/mol")
    _print(f"  hairpin_dg_buffer         : {hairpin_dg_buffer} kcal/mol")
    _print(f"  homodimer_dg_threshold    : {homodimer_dg_threshold} kcal/mol")
    _print(f"  homodimer_dg_buffer       : {homodimer_dg_buffer} kcal/mol")
    _print(f"  gc_clamp_length           : {gc_clamp_length} bp")
    _print(f"  gc_clamp_min_gc           : {gc_clamp_min_gc}")
    _print(f"  heterodimer_dg_threshold  : {heterodimer_dg_threshold} kcal/mol")
    _print(f"  heterodimer_dg_buffer     : {heterodimer_dg_buffer} kcal/mol")
    _print("  -- Misprime / heterodimer detection -------------------------")
    _print(f"  misprime_length           : {misprime_length}")
    _print(f"  misprime_tip_length       : {misprime_tip_length}")
    _print(f"  misprime_max_mismatches   : {misprime_max_mismatches}")
    _print(f"  misprime_kmer_size        : {misprime_kmer_size}")
    _print(f"  hetdimer_kmer_size        : {hetdimer_kmer_size}")
    _print(f"  heterodimer_neighbor_range: {heterodimer_neighbor_range}")
    _print("  -- DP penalty weights ----------------------------------------")
    _print(f"  penalty_tm                : {penalty_tm}")
    _print(f"  penalty_gc                : {penalty_gc}")
    _print(f"  penalty_homopolymer       : {penalty_homopolymer}")
    _print(f"  penalty_dinucleotide      : {penalty_dinucleotide}")
    _print(f"  penalty_hairpin           : {penalty_hairpin}")
    _print(f"  penalty_homodimer         : {penalty_homodimer}")
    _print(f"  penalty_gc_clamp          : {penalty_gc_clamp}")
    _print(f"  penalty_misprime          : {penalty_misprime}")
    _print(f"  penalty_hetdimer_map      : {penalty_hetdimer_map}")
    _print(f"  penalty_body_homopolymer  : {penalty_body_homopolymer}")
    _print(f"  penalty_body_dinucleotide : {penalty_body_dinucleotide}")
    _print(f"  penalty_severity_mult     : {penalty_severity_multiplier}")
    _print("  -- Fragment cut scoring --------------------------------------")
    _print(f"  cut_search_window         : {cut_search_window}")
    _print(f"  cut_score_window          : {cut_score_window}")
    _print(f"  cut_score_repeat_penalty  : {cut_score_repeat_penalty}")
    _print("  -- Junction and outer primers --------------------------------")
    _print(f"  junction_overlap_length   : {junction_overlap_length}")
    _print(f"  junction_tm_range         : {junction_tm_range}")
    _print(f"  outer_primer_length_range : {outer_primer_length_range}")
    _print(f"  outer_primer_tm_range     : {outer_primer_tm_range}")
    _print("  -- Assembly method -------------------------------------------")
    _print(f"  assembly_method           : {assembly_method}")
    _print(f"  goldengate_enzyme         : {goldengate_enzyme}")
    _print("  -- Thermodynamic conditions -----------------------------------")
    _print(f"  mv_conc                   : {mv_conc} mM Na+")
    _print(f"  dv_conc                   : {dv_conc} mM Mg2+")
    _print(f"  dntp_conc                 : {dntp_conc} mM")
    _print(f"  dna_conc                  : {dna_conc} nM")
    _print("  -- Workflow control ------------------------------------------")
    _print(f"  run_oligos                : {run_oligos}")
    _print(f"  run_primers               : {run_primers}")
    _print(f"  fragment_indices          : {fragment_indices}")
    _print(f"  n_workers                 : {n_workers}")
    _print(f"  output_csv                : {output_csv}")
    _print(f"  output_log                : {output_log}")
    _print("=" * 60)
    _print()

    # Input cleaning

    # Supported Formats: raw ACGT, single-entry FASTA, and multi-line GenBank format
    sequence = clean_input_sequence(sequence)

    # Basic input validation
    invalid = set(sequence) - set('ACGT')
    if invalid:
        _raise(ValueError(
            f"Sequence contains invalid characters: {invalid}, Only A, T, G, C are supported (no IUPAC ambiguity codes)"
        ))
    if len(sequence) < oligo_length_range[0]:
        _raise(ValueError(
            f"Sequence length ({len(sequence)} bp) is shorter than the minimum oligo length ({oligo_length_range[0]} bp), Nothing to tile"
        ))


    # Cross-parameter logic checks
    min_ol, max_ol = oligo_length_range
    min_ov, max_ov = overlap_length_range

    if min_ol < 1:
        _raise(ValueError(f"oligo_length_range min ({min_ol}) must be >= 1 bp."))
    if min_ol > max_ol:
        _raise(ValueError(f"oligo_length_range: min ({min_ol}) must be <= max ({max_ol})."))
    if min_ov < 1:
        _raise(ValueError(f"overlap_length_range min ({min_ov}) must be >= 1 bp."))
    if min_ov > max_ov:
        _raise(ValueError(f"overlap_length_range: min ({min_ov}) must be <= max ({max_ov})."))
    if min_ov >= min_ol:
        _raise(ValueError(
            f"Impossible: minimum overlap ({min_ov} bp) must be strictly less than minimum oligo length ({min_ol} bp)"
        ))
    if overlap_tm_range[0] > overlap_tm_range[1]:
        _raise(ValueError(
            f"overlap_tm_range: min ({overlap_tm_range[0]}C) must be <= max ({overlap_tm_range[1]}C)"
        ))
    if overlap_gc_range[0] > overlap_gc_range[1]:
        _raise(ValueError(
            f"overlap_gc_range: min ({overlap_gc_range[0]}%) must be <= max ({overlap_gc_range[1]}%)"
        ))
    if fragment_length_range:
        if fragment_length_range[0] > fragment_length_range[1]:
            _raise(ValueError("fragment_length_range: min must be <= max."))
        if fragment_length_range[0] <= max_ov:
            _raise(ValueError(
                f"Impossible: minimum fragment length ({fragment_length_range[0]} bp) must be greater than maximum overlap ({max_ov} bp)"
            ))
        if max_ol > fragment_length_range[1]:
            _raise(ValueError(
                f"oligo_length_range max ({max_ol} bp) exceeds fragment_length_range max ({fragment_length_range[1]} bp). An oligo cannot be longer than the fragment it tiles."
            ))
        if min_ol > fragment_length_range[1]:
            _raise(ValueError(
                f"oligo_length_range min ({min_ol} bp) exceeds fragment_length_range max ({fragment_length_range[1]} bp). Even the shortest oligo cannot fit in the longest fragment."
            ))

    # Overlap can never be the entire oligo, at least 1 non-overlap base required
    if max_ov >= min_ol:
        _raise(ValueError(
            f"overlap_length_range max ({max_ov} bp) must be strictly less than oligo_length_range min ({min_ol} bp). Every oligo must have at least 1 non-overlap base."
        ))

    # Warn if overlap_length_range contains sub-10 bp values
    if min_ov < MIN_RELIABLE_TM_LENGTH:
        _warnings.warn(
            f"overlap_length_range minimum ({min_ov} bp) is below {MIN_RELIABLE_TM_LENGTH} bp. The SantaLucia Tm model may be unreliable for such short overlaps. Consider raising overlap_length_range[0] to at least {MIN_RELIABLE_TM_LENGTH}",
            UserWarning,
            stacklevel=2,
        )

    # Assembly method validation and enforcement
    _am = assembly_method.lower().strip() # normalise string
    if _am not in ASSEMBLY_METHODS:
        _raise(ValueError(
            f"assembly_method must be one of {ASSEMBLY_METHODS}, Got: '{assembly_method}'"
        ))

    if _am == "golden_gate":
        if junction_overlap_length != 0:
            _warnings.warn(
                f"Golden Gate Assembly selected: junction_overlap_length ({junction_overlap_length} bp) is ignored. Golden Gate uses 4-bp sticky ends with zero sequence overlap between fragments. junction_overlap_length has been set to 0.",
                UserWarning,
                stacklevel=2,
            )
        # GG uses 4-bp sticky ends: force junction_overlap to 0
        junction_overlap_length = 0

        # naively estimate fragment count, exit early if it exceeds the GG overhang pool
        if fragment_length_range is not None:
            _midpoint      = (fragment_length_range[0] + fragment_length_range[1]) / 2
            _est_frags     = max(1, round(len(sequence) / _midpoint))
            if _est_frags > MAX_GOLDENGATE_FRAGMENTS:
                _raise(ValueError(
                    f"Golden Gate infeasible: estimated fragment count ({_est_frags}) exceeds the maximum of {MAX_GOLDENGATE_FRAGMENTS} fragments supported by the orthogonal overhang pool. Consider splitting the construct into sub-assemblies or using Gibson Assembly instead"
                ))

        # Enzyme Restriction site scan
        _enzyme_order   = [goldengate_enzyme] + [
            e for e in GOLDEN_GATE_ENZYMES if e != goldengate_enzyme
        ]
        _gg_scan_results = {}
        _active_enzyme   = None
        for _enz in _enzyme_order:
            _hits = _scan_restriction_sites(sequence, _enz)
            _gg_scan_results[_enz] = _hits
            if len(_hits) == 0:
                _active_enzyme = _enz
                break

        _print(f"  Restriction site scan (full input sequence):")
        for _enz, _hits in _gg_scan_results.items():
            if _enz == _active_enzyme:
                _status = "PASS (selected)"
            else:
                _status = f"FAIL ({len(_hits)} site(s) at positions {_hits})"
            _print(f"    {_enz:6s}: {_status}")

        if _active_enzyme is None:
            _bsai  = _gg_scan_results.get("BsaI",  [])
            _bsmbi = _gg_scan_results.get("BsmBI", [])
            _bbsi  = _gg_scan_results.get("BbsI",  [])
            _raise(ValueError(
                f"Golden Gate infeasible: all three Type IIS enzymes (BsaI, BsmBI, BbsI) have internal recognition sites \n"
                f"BsaI  : {len(_bsai)} site(s) at positions {_bsai}\n"
                f"BsmBI : {len(_bsmbi)} site(s) at positions {_bsmbi}\n"
                f"BbsI  : {len(_bbsi)} site(s) at positions {_bbsi}\n"
                f"Options include using Gibson Assembly or silently mutating recognition sites via codon optimisation"   
            ))

        if _active_enzyme != goldengate_enzyme:
            _warnings.warn(
                f"Golden Gate: chosen enzyme '{goldengate_enzyme}' has internal recognition sites in the input sequence. Automatically switched to '{_active_enzyme}' which has no internal sites",
                UserWarning,
                stacklevel=2,
            )
        goldengate_enzyme = _active_enzyme


    if _am == "gibson":
        # raise junction overlap to 40bp minimum (T5 exonuclease chew-back safety)
        if junction_overlap_length < GIBSON_MIN_JUNCTION_OVERLAP:
            _warnings.warn(
                f"Gibson Assembly: junction_overlap_length ({junction_overlap_length} bp) raised to {GIBSON_MIN_JUNCTION_OVERLAP} bp (T5 exonuclease minimum)",
                UserWarning,
                stacklevel=2,
            )
            junction_overlap_length = GIBSON_MIN_JUNCTION_OVERLAP

        # relax junction Tm minimum to 48°C (Gibson runs at 50°C, PCR default 58°C flags valid overlaps)
        if junction_tm_range == DEFAULT_JUNCTION_TM_RANGE:
            junction_tm_range = (GIBSON_MIN_JUNCTION_TM, junction_tm_range[1])
            _warnings.warn(
                f"Gibson Assembly: junction Tm minimum relaxed from {DEFAULT_JUNCTION_TM_RANGE[0]}°C to {GIBSON_MIN_JUNCTION_TM}°C (Gibson runs at 50°C, PCR default is too strict)",
                UserWarning,
                stacklevel=2,
            )

    # Penalties dict
    penalties = {
        'tm':                  penalty_tm,
        'gc':                  penalty_gc,
        'homopolymer':         penalty_homopolymer,
        'dinucleotide':        penalty_dinucleotide,
        'hairpin':             penalty_hairpin,
        'homodimer':           penalty_homodimer,
        'gc_clamp':            penalty_gc_clamp,
        'misprime':            penalty_misprime,
        'hetdimer_map':        penalty_hetdimer_map,
        'body_homopolymer':    penalty_body_homopolymer,
        'body_dinucleotide':   penalty_body_dinucleotide,
        'severity_multiplier': penalty_severity_multiplier,
    }

    # Startup header
    seq_label = output_csv.rsplit('.', 1)[0] if output_csv else 'unnamed'
    _print("\n" + "=" * 60)
    _print("  OligoFit")
    _print("=" * 60)
    _print(f"  Sequence              : {seq_label}  ({len(sequence)} bp)")
    _print(f"  Fragment length       : {fragment_length_range[0]}–{fragment_length_range[1]} bp"
          if fragment_length_range else f"  Fragment length       : single fragment (no split)")
    _print(f"  Oligo length          : {oligo_length_range[0]}–{oligo_length_range[1]} bp")
    _print(f"  Overlap length        : {overlap_length_range[0]}–{overlap_length_range[1]} bp")
    _print(f"  Overlap Tm target     : {overlap_tm_range[0]:.0f}–{overlap_tm_range[1]:.0f} °C")
    _print(f"  --- Key penalty weights ---")
    _print(f"  Tm                    : {penalty_tm}   "
          f"Hairpin dG threshold  : {hairpin_dg_threshold} kcal/mol")
    _print(f"  Homodimer             : {penalty_homodimer}   "
          f"Homodimer dG threshold: {homodimer_dg_threshold} kcal/mol")
    _print(f"  Misprime              : {penalty_misprime}   "
          f"Heterodimer dG thresh : {heterodimer_dg_threshold} kcal/mol")
    _print(f"  Hetdimer map          : {penalty_hetdimer_map}   "
          f"Severity multiplier   : {penalty_severity_multiplier}")
    _print("=" * 60)

    _am_label = {
        'stitching_pcr': 'Stitching PCR',
        'gibson':        'Gibson Assembly',
        'golden_gate':   'Golden Gate',
    }.get(_am, _am)
    _print(f"  Assembly method       : {_am_label}")
    if _am == "gibson":
        _print(f"  Junction overlap      : {junction_overlap_length} bp  "
              f"(Gibson min: {GIBSON_MIN_JUNCTION_OVERLAP} bp)")
        _print(f"  Junction Tm min       : {junction_tm_range[0]}°C  "
              f"(Gibson isothermal reaction: 50°C)")
    _print("=" * 60)

    # --- Sequence complexity assessment (runs before any tiling) ---
    _seq_assessment = evaluate_dna_sequence(
        seq=sequence,
        max_homopolymer=max_homopolymer_run,
        max_dinucleotide=max_dinucleotide_repeat,
    )
    _print("\n  -- Sequence Complexity Assessment --")
    _print(f"  Length                : {len(sequence)} bp")
    _print(f"  GC content            : {_seq_assessment['gc_percent']:.1f}%")
    _hp = _seq_assessment['homopolymers']
    _dn = _seq_assessment['dinucleotide_repeats']
    _print(f"  Homopolymer runs (>={max_homopolymer_run}bp) : {_hp['total_count']}"
           + (f"  {_hp['details']}" if _hp['details'] else "  none"))
    _print(f"  Dinucleotide repeats  : {_dn['total_count']}"
           + (f"  {_dn['details']}" if _dn['details'] else "  none"))
    _print("=" * 60)

    _print(f"\n[1/6] Fragmenting sequence …")
    gg_internal_overhangs = []
    if fragment_length_range is None:
        fragments = [sequence]
    elif _am == "golden_gate":
        fragments, gg_internal_overhangs = _fragment_sequence_golden_gate(
            sequence=sequence,
            frag_range=fragment_length_range,
            cut_search_window=cut_search_window,
            max_homopolymer=max_homopolymer_run,
            max_dinucleotide=max_dinucleotide_repeat,
            cut_score_window=cut_score_window,
            cut_score_repeat_penalty=cut_score_repeat_penalty,
            enzyme_name=goldengate_enzyme,
        )
    else:
        fragments = _fragment_sequence(
            sequence, fragment_length_range,
            junction_overlap=junction_overlap_length,
            cut_search_window=cut_search_window,
            max_homopolymer=max_homopolymer_run,
            max_dinucleotide=max_dinucleotide_repeat,
            cut_score_window=cut_score_window,
            cut_score_repeat_penalty=cut_score_repeat_penalty,
        )
    _print(f"      {len(fragments)} fragment(s): {[len(f) for f in fragments]} bp")

    # fragment_indices filter
    if fragment_indices is not None:
        valid = set(range(1, len(fragments) + 1))
        requested = set(fragment_indices)
        bad = requested - valid
        if bad:
            _warnings.warn(
                f"fragment_indices contains out-of-range values {sorted(bad)}. "
                f"Valid range is 1–{len(fragments)}. These will be ignored.",
                UserWarning, stacklevel=2,
            )
        active_indices = sorted(requested & valid)
        if not active_indices:
            _raise(ValueError(
                "fragment_indices produced no valid fragments to tile. "
                f"Valid range is 1–{len(fragments)}."
            ))
    else:
        active_indices = list(range(1, len(fragments) + 1))  # all, 1-based

    # Step 1b: Build risk maps
    _print(f"[2/6] Building misprime (k={misprime_kmer_size}) and heterodimer RC (k={hetdimer_kmer_size}) risk maps …")
    frag_risk_maps = [
        _build_misprime_risk_map(sequence=frag, kmer_size=misprime_kmer_size)
        for frag in fragments
    ]

    frag_hetdimer_maps = [
        _build_hetdimer_risk_map(fragment=frag, kmer_size=hetdimer_kmer_size)
        for frag in fragments
    ]

    # Global start positions: each non-last fragment advances len(frag) - junction_overlap bp
    frag_global_starts = []
    cumulative = 0
    for i, frag in enumerate(fragments):
        frag_global_starts.append(cumulative)
        if i < len(fragments) - 1:
            cumulative += len(frag) - junction_overlap_length
        else:
            cumulative += len(frag)

    # Step 2: Tile fragments
    _gg_overhangs = None
    if _am == "golden_gate":
        if len(fragments) > MAX_GOLDENGATE_FRAGMENTS:
            _raise(ValueError(
                f"Golden Gate infeasible: {len(fragments)} fragments produced, "
                f"exceeding the maximum of {MAX_GOLDENGATE_FRAGMENTS} fragments."
            ))
        _gg_overhangs = [None] + gg_internal_overhangs + [None]
        _print(f"      Native internal overhangs chosen: {len(gg_internal_overhangs)}  "
              f"(for {len(fragments)} fragments)")

    if not run_oligos:
        _print("[3/6] Oligo tiling skipped (run_oligos=False).")
        _print("[4/6] Post-processing skipped.")
        frag_df = pd.DataFrame({
            'fragment_num':       list(range(1, len(fragments) + 1)),
            'fragment_length_bp': [len(f) for f in fragments],
            'fragment_sequence':  fragments,
        })
        _print("\n  Fragment table returned. No oligo, junction, or primer design run.")
        _write_log()
        return pd.DataFrame(), frag_df, pd.DataFrame()

    n_active = len(active_indices)
    _print(f"[3/6] Tiling {n_active} fragment(s) using A* DP …")

    worker_args = []
    
    gg_rec, gg_spacer = "", ""
    if _am == "golden_gate":
        gg_rec    = GOLDEN_GATE_ENZYMES[goldengate_enzyme]["recognition"]
        gg_spacer = GOLDEN_GATE_ENZYMES[goldengate_enzyme]["fwd_spacer"]

    for frag_idx, fragment in enumerate(fragments):
        if (frag_idx + 1) not in active_indices:
            continue
            
        fwd_cassette, rev_cassette = "", ""
        if _am == "golden_gate":
            left_ov = _gg_overhangs[frag_idx]
            right_ov = _gg_overhangs[frag_idx + 1]
            if left_ov is not None:
                fwd_cassette = GOLDEN_GATE_EXTRA_BASES + gg_rec + gg_spacer
            if right_ov is not None:
                # Reverse cassette appended to the right of the sequence
                # For top strand, it is the reverse complement of the forward cassette format
                rev_cassette = reverse_complement(GOLDEN_GATE_EXTRA_BASES + gg_rec + gg_spacer)
                
        _print(f"      Fragment {frag_idx+1}/{len(fragments)}, {len(fragment)} bp, starting …")
        worker_args.append((
            fragment, frag_idx,
            oligo_length_range, overlap_length_range,
            overlap_tm_range, tm_warning_buffer, overlap_gc_range,
            max_homopolymer_run, max_dinucleotide_repeat,
            repeat_vicinity_buffer, hairpin_dg_threshold,
            mv_conc, dv_conc, dntp_conc, dna_conc,
            fragment_length_range,
            frag_risk_maps[frag_idx], frag_global_starts[frag_idx],
            penalties, homodimer_dg_threshold,
            gc_clamp_length, gc_clamp_min_gc,
            misprime_kmer_size,
            frag_hetdimer_maps[frag_idx], hetdimer_kmer_size,
            gc_warning_buffer, hairpin_dg_buffer, homodimer_dg_buffer,
            fwd_cassette, rev_cassette
        ))

    frag_results: dict = {}

    try:
        effective_workers = 1 if len(worker_args) == 1 else n_workers
        if effective_workers == 1:
            raise ValueError("single worker, use sequential path")

        completed = 0
        with ProcessPoolExecutor(max_workers=effective_workers) as executor:
            futures = {
                executor.submit(_tile_fragment_worker, args): args[1]
                for args in worker_args
            }
            for future in as_completed(futures):
                frag_idx, records = future.result()
                frag_results[frag_idx] = records
                completed += 1
                _print(f"      Fragment {frag_idx+1}/{len(fragments)}, done "
                      f"({len(records)} oligos)  [{completed}/{n_active} complete]")

    except (concurrent.futures.BrokenExecutor, OSError, ValueError):
        for args in worker_args:
            frag_idx, records = _tile_fragment_worker(args)
            frag_results[frag_idx] = records
            _print(f"      Fragment {frag_idx+1}/{len(fragments)}, done "
                  f"({len(records)} oligos)  [{len(frag_results)}/{n_active} complete]")

    # Reconstruct in fragment order and assign global_oligo_num + global coords
    all_records = []
    global_num  = 1
    for frag_idx in range(len(fragments)):
        if (frag_idx + 1) not in set(active_indices):
            continue
            
        _validate_fragment_assembly(fragments[frag_idx], frag_results[frag_idx])
        _print(f"      [OK] Fragment {frag_idx+1} properly re-assembles from {len(frag_results[frag_idx])} oligos")
        
        frag_offset = frag_global_starts[frag_idx]
        for rec in frag_results[frag_idx]:
            rec['global_oligo_num'] = global_num
            rec['global_start'] = frag_offset + rec['frag_start']
            rec['global_end']   = frag_offset + rec['frag_end']
            global_num += 1
            all_records.append(rec)

    # For Golden Gate, build the cassettes directly into terminal oligos.
    if _am == "golden_gate" and _gg_overhangs is not None:
        _print("      Embedding the cassettes into terminal oligos (Golden Gate) ...")
        enzyme  = GOLDEN_GATE_ENZYMES[goldengate_enzyme]
        rec_seq = enzyme["recognition"]     
        rc_rec  = enzyme["rc_recognition"]  
        spacer  = enzyme["fwd_spacer"]      
        extras  = GOLDEN_GATE_EXTRA_BASES   

        # Group all_records indices by fragment number
        frag_rec_indices: dict = {}
        for i, rec in enumerate(all_records):
            fnum = rec['fragment']
            frag_rec_indices.setdefault(fnum, []).append(i)

        for fnum in sorted(frag_rec_indices):
            idxs     = frag_rec_indices[fnum]
            frag_idx = fnum - 1
            left_ov  = _gg_overhangs[frag_idx]       # 4-bp left overhang
            right_ov = _gg_overhangs[frag_idx + 1]   # 4-bp right overhang

            # First oligo
            # Prepend: [GCATAC] - [rec] - [spacer]
            if left_ov is not None:
                first = all_records[idxs[0]]
                left_cassette = extras + rec_seq + spacer
                first['oligo_sequence'] = left_cassette + first['oligo_sequence']
                first['oligo_length']   = len(first['oligo_sequence'])

            # Last oligo
            # Top strand architecture at the right end: [fragment] - [spacer] - [rc_rec] - [GCATAC]
            if right_ov is not None:
                last = all_records[idxs[-1]]
                right_cassette_top = spacer + rc_rec + extras
                if last['direction'] == 'forward':
                    last['oligo_sequence'] = last['oligo_sequence'] + right_cassette_top
                else:
                    last['oligo_sequence'] = reverse_complement(right_cassette_top) + last['oligo_sequence']
                last['oligo_length'] = len(last['oligo_sequence'])

    # Step 2b: Misprime check
    _print(f"[4/6] Running post-tiling misprime scan ({len(all_records)} oligos) …")
    misprime_flags, misprime_warns = check_misprime_all(
        records=all_records,
        fragments=fragments,
        misprime_length=misprime_length,
        tip_length=misprime_tip_length,
        max_mismatches=misprime_max_mismatches,
        frag_global_starts=frag_global_starts,
    )

    for i, rec in enumerate(all_records):
        rec['flag_misprime'] = misprime_flags[i]
        if misprime_warns[i]:
            existing = rec.get('warnings', 'OK')
            rec['warnings'] = (misprime_warns[i] if existing == 'OK'
                               else existing + '; ' + misprime_warns[i])

    # Step 2c: Heterodimer check
    _print(f"      Running heterodimer check (±{heterodimer_neighbor_range} neighbors) …")
    hetdimer_flags, hetdimer_warns = check_heterodimers_all(
        records=all_records,
        dg_threshold=heterodimer_dg_threshold,
        neighbor_range=heterodimer_neighbor_range,
        dg_buffer=heterodimer_dg_buffer,
        mv_conc=mv_conc, dv_conc=dv_conc,
        dntp_conc=dntp_conc, dna_conc=dna_conc,
    )

    for i, rec in enumerate(all_records):
        rec['flag_heterodimer'] = hetdimer_flags[i]
        if hetdimer_warns[i]:
            existing = rec.get('warnings', 'OK')
            rec['warnings'] = (hetdimer_warns[i] if existing == 'OK'
                               else existing + '; ' + hetdimer_warns[i])

    # Step 3: Assemble the DataFrame
    _print("[5/6] Assembling output DataFrame …")
    col_order = [
        'global_oligo_num', 'oligo_name', 'fragment', 'oligo_in_frag',
        'direction', 'frag_start', 'frag_end', 'global_start', 'global_end',
        'oligo_length', 'oligo_gc_%',
        'oligo_sequence', 'overlap_sequence', 'overlap_length',
        'overlap_tm_C', 'overlap_gc_%',
        'flag_tm_out_of_range', 'flag_overlap_gc', 'flag_homopolymer',
        'flag_dinucleotide_repeat', 'flag_hairpin', 'flag_homodimer',
        'flag_gc_clamp', 'flag_heterodimer', 'flag_misprime',
        'flag_body_homopolymer', 'flag_body_dinucleotide',
        'penalty_score', 'is_last_in_frag', 'warnings',
    ]
    df = pd.DataFrame(all_records)[col_order]

    # Step 4: Design fragment table (includes junction info)
    _print("[6/6] Designing fragment table and outer primers …")
    fragment_df = _design_fragment_table(
        fragments=fragments,
        junction_overlap=junction_overlap_length,
         junction_tm_range=junction_tm_range,
        tm_warning_buffer=tm_warning_buffer,
        gc_range=overlap_gc_range,
        max_homopolymer=max_homopolymer_run,
        max_dinucleotide=max_dinucleotide_repeat,
        hairpin_dg_threshold=hairpin_dg_threshold,
        mv_conc=mv_conc, dv_conc=dv_conc,
        dntp_conc=dntp_conc, dna_conc=dna_conc,
        frag_global_starts=frag_global_starts,
        assembly_method=_am,
    )

    if _am == "golden_gate" and _gg_overhangs is not None:
        # Update fragment table sequences to include embedded Golden Gate cassettes
        enzyme  = GOLDEN_GATE_ENZYMES[goldengate_enzyme]
        rec_seq = enzyme["recognition"]     
        rc_rec  = enzyme["rc_recognition"]  
        spacer  = enzyme["fwd_spacer"]      
        extras  = GOLDEN_GATE_EXTRA_BASES   

        new_frag_seqs = []
        new_frag_lens = []
        for i, row in fragment_df.iterrows():
            frag_idx = row['fragment_num'] - 1
            raw_seq  = row['fragment_sequence']
            
            left_ov  = _gg_overhangs[frag_idx]
            right_ov = _gg_overhangs[frag_idx + 1]
            
            left_cassette = ""
            if left_ov is not None:
                left_cassette = extras + rec_seq + spacer
                
            right_cassette = ""
            if right_ov is not None:
                right_cassette = spacer + rc_rec + extras
            
            full_seq = left_cassette + raw_seq + right_cassette
            new_frag_seqs.append(full_seq)
            new_frag_lens.append(len(full_seq))
            
        fragment_df['fragment_sequence'] = new_frag_seqs
        fragment_df['fragment_length_bp'] = new_frag_lens


    # Step 4b: Outer PCR primers
    if run_primers:
        _base_outer_df = _design_outer_primers(

            fragments=fragment_df['fragment_sequence'].tolist(),
            primer_length_range=outer_primer_length_range,
            primer_tm_range=outer_primer_tm_range,
            tm_warning_buffer=tm_warning_buffer,
            gc_range=overlap_gc_range,  # overlap_gc_range reused here
            max_homopolymer=max_homopolymer_run,
            hairpin_dg_threshold=hairpin_dg_threshold,
            mv_conc=mv_conc, dv_conc=dv_conc,
            dntp_conc=dntp_conc, dna_conc=dna_conc,
        )
        outer_df = _base_outer_df
    else:
        _print("      Outer primer design skipped (run_primers=False).")
        outer_df = pd.DataFrame()

    _validate_full_assembly(
        original_seq=sequence,
        fragment_df=fragment_df,
        assembly_method=_am,
        junction_overlap_length=junction_overlap_length,
        gg_overhangs=_gg_overhangs,
        goldengate_enzyme=goldengate_enzyme
    )
    _print(f"      [OK] Original Sequence properly re-assembled from {len(fragment_df)} fragments")

    # Step 5: Print summary
    n_warnings       = (df['warnings'] != 'OK').sum()
    valid_tms        = df['overlap_tm_C'].dropna()
    penalty_scores   = df['penalty_score']
    n_gc_flags          = df['flag_overlap_gc'].sum()
    n_hp_flags          = df['flag_homopolymer'].sum()
    n_dn_flags          = df['flag_dinucleotide_repeat'].sum()
    n_hairpin_flags     = df['flag_hairpin'].sum()
    n_homodimer_flags   = df['flag_homodimer'].sum()
    n_gc_clamp_flags    = df['flag_gc_clamp'].sum()
    n_hetdimer_flags    = df['flag_heterodimer'].sum()
    n_tm_flags          = df['flag_tm_out_of_range'].sum()
    n_misprime_flags    = df['flag_misprime'].sum()
    n_body_hp_flags     = df['flag_body_homopolymer'].sum()
    n_body_dn_flags     = df['flag_body_dinucleotide'].sum()

    _print("\n" + "=" * 60)
    _print("  OligoFit SUMMARY")
    _print("=" * 60)
    _print(f"  Input sequence        : {seq_label}  ({len(sequence)} bp)")
    _print(f"  Fragments             : {len(fragments)}  "
          f"({fragment_df['fragment_length_bp'].tolist() if not fragment_df.empty else [len(f) for f in fragments]} bp)")
    _print(f"  Algorithm             : A* DP (globally optimal)")
    _print(f"  Risk maps             : misprime k={misprime_kmer_size}, "
          f"heterodimer RC k={hetdimer_kmer_size}  (fragment-local)")
    if len(fragments) > 1 and junction_overlap_length > 0:
        _print(f"  Junction overlap      : {junction_overlap_length} bp "
              f"(shared between adjacent fragments for stitching PCR)")
    _print(f"  Total oligos          : {len(df)}")
    if len(valid_tms):
        _print(f"  Overlap Tm range      : {valid_tms.min():.1f}C -- "
              f"{valid_tms.max():.1f}C")
        _print(f"  Overlap Tm mean       : {valid_tms.mean():.1f}C")
        _print(f"  Overlap Tm std dev    : {valid_tms.std():.2f}C "
              f"")
    _print(f"  --- Assembly junction quality flags ---")
    _print(f"  Tm out of range       : {n_tm_flags} oligo(s)")
    _print(f"  Overlap GC out range  : {n_gc_flags} oligo(s)")
    _print(f"  Homopolymer (overlap) : {n_hp_flags} oligo(s)")
    _print(f"  Dinucleotide (overlap): {n_dn_flags} oligo(s)")
    _print(f"  Potential hairpins    : {n_hairpin_flags} oligo(s)")
    _print(f"  Homodimer risk        : {n_homodimer_flags} oligo(s)")
    _print(f"  Weak GC clamp         : {n_gc_clamp_flags} oligo(s)")
    _print(f"  Heterodimer risk      : {n_hetdimer_flags} oligo(s)")
    _print(f"  Potential misprimes   : {n_misprime_flags} oligo(s)")
    _print(f"  --- Synthesis quality flags (oligo body) ---")
    _print(f"  Homopolymer (body)    : {n_body_hp_flags} oligo(s)  [synthesis / template-slippage risk]")
    _print(f"  Dinucleotide (body)   : {n_body_dn_flags} oligo(s)  [synthesis / template-slippage risk]")
    _print(f"  Total with warnings   : {n_warnings} oligo(s)")
    _print(f"  Mean penalty score    : {penalty_scores.mean():.3f}  "
          f"(max: {penalty_scores.max():.3f}, worst: "
          f"{df.loc[penalty_scores.idxmax(), 'oligo_name'] if len(df) else 'N/A'})")
    if n_warnings == 0:
        _print(f"  Status                : OK, all oligos passed all checks")
    else:
        _print(f"  Status                : See 'warnings' column for details")

    # Step 5b: NOTES section -- actionable recommendations based on flag analysis
    n_total_oligos = len(df)
    if n_total_oligos > 0:
        # Count oligos that have at least one flag (any flag column is True)
        _flag_cols = [
            'flag_tm_out_of_range', 'flag_overlap_gc', 'flag_homopolymer',
            'flag_dinucleotide_repeat', 'flag_hairpin', 'flag_homodimer',
            'flag_gc_clamp', 'flag_heterodimer', 'flag_misprime',
            'flag_body_homopolymer', 'flag_body_dinucleotide',
        ]
        _existing_flag_cols = [c for c in _flag_cols if c in df.columns]
        _n_flagged_oligos = int((df[_existing_flag_cols].any(axis=1)).sum())
        _flagged_ratio = _n_flagged_oligos / n_total_oligos

        # Build per-flag count map (only flags that fired at least once)
        _flag_counts = {
            col: int(df[col].sum())
            for col in _existing_flag_cols
            if int(df[col].sum()) > 0
        }

        # Human-readable flag names
        _flag_display_names = {
            'flag_tm_out_of_range':     'Tm out of range',
            'flag_overlap_gc':          'Overlap GC out of range',
            'flag_homopolymer':         'Homopolymer (overlap)',
            'flag_dinucleotide_repeat': 'Dinucleotide repeat (overlap)',
            'flag_hairpin':             'Hairpin',
            'flag_homodimer':           'Homodimer',
            'flag_gc_clamp':            'Weak GC clamp',
            'flag_heterodimer':         'Heterodimer',
            'flag_misprime':            'Misprime',
            'flag_body_homopolymer':    'Homopolymer (body)',
            'flag_body_dinucleotide':   'Dinucleotide repeat (body)',
        }

        _has_notes = False

        if _flagged_ratio >= 0.80 or _flag_counts:
            _print(f"\n  --- NOTES ---")
            _has_notes = True

        # Note 1: High flag ratio warning
        if _flagged_ratio >= 0.80:
            _pct = round(_flagged_ratio * 100, 1)
            _print(f"  [!] {_pct}% of oligos ({_n_flagged_oligos}/{n_total_oligos}) carry "
                   f"at least one flag.")
            _print(f"      Consider expanding the oligo length range (currently "
                   f"{oligo_length_range[0]}-{oligo_length_range[1]} bp) and/or the")
            _print(f"      overlap length range (currently "
                   f"{overlap_length_range[0]}-{overlap_length_range[1]} bp) so the "
                   f"algorithm can")
            _print(f"      explore more placement options and find lower-penalty oligo "
                   f"configurations.")

        # Note 2: Identify the highest (most common) flag type(s)
        if _flag_counts:
            _max_count = max(_flag_counts.values())
            _top_flags = [col for col, cnt in _flag_counts.items() if cnt == _max_count]
            _top_names = [_flag_display_names.get(f, f) for f in _top_flags]

            if len(_top_names) == 1:
                _print(f"  [i] The most prevalent flag is '{_top_names[0]}' "
                       f"({_max_count} oligo(s)).")
            else:
                _joined = "' and '".join(_top_names)
                _print(f"  [i] The most prevalent flags are '{_joined}' "
                       f"({_max_count} oligo(s) each).")
            _print(f"      Consider raising its penalty weight and/or expanding the oligo "
                   f"and overlap")
            _print(f"      length ranges so the algorithm has a wider search space to find")
            _print(f"      placements that avoid this issue.")

    if not fragment_df.empty:
        # Rows with a junction (all except the last fragment)
        junc_rows = fragment_df[fragment_df['junction_to_frag'].notna()]
        n_junc_warnings = (junc_rows['warnings'] != 'OK').sum() if not junc_rows.empty else 0
        _print(f"\n  --- Fragment table ---")
        _print(f"  Fragments             : {len(fragment_df)}")
        if not junc_rows.empty:
            _print(f"  Junctions             : {len(junc_rows)}")
            for _, frow in junc_rows.iterrows():
                status = 'OK' if frow['warnings'] == 'OK' else 'FLAGGED'
                tm_s   = f"{frow['junction_tm_C']:.1f}C" if frow['junction_tm_C'] is not None else 'N/A'
                gc_s   = f"{frow['junction_gc_%']:.1f}%" if frow['junction_gc_%'] is not None else 'N/A'
                label  = f"Frag{int(frow['fragment_num']):02d}→Frag{int(frow['junction_to_frag']):02d}"
                _print(f"  {label:20s}  Tm={tm_s}  GC={gc_s}  [{status}]")
            if n_junc_warnings > 0:
                _print(f"  Junction warnings     : {n_junc_warnings}, "
                      f"see fragment CSV for details")
            else:
                _print(f"  Junction status       : OK, all junctions passed all checks")

    if not outer_df.empty:
        n_outer_warn = (outer_df['warnings'] != 'OK').sum()
        _print(f"\n  --- Outer PCR primer summary (via primer3.design_primers) ---")
        _print(f"  Primers designed      : {len(outer_df)} "
              f"({len(fragments)} fragment(s) × 2)")
        valid_outer_tms = outer_df['primer_tm_C'].dropna()
        if len(valid_outer_tms):
            _print(f"  Primer Tm range       : {valid_outer_tms.min():.1f}C -- "
                  f"{valid_outer_tms.max():.1f}C")
        for _, prow in outer_df.iterrows():
            status  = 'OK' if prow['warnings'] == 'OK' else 'FLAGGED'
            tm_str  = f"{prow['primer_tm_C']:.1f}C" if prow['primer_tm_C'] is not None else 'N/A'
            gc_str  = f"{prow['primer_gc_%']:.1f}%" if prow['primer_gc_%'] is not None else 'N/A'
            len_str = f"{prow['primer_length']}bp" if prow['primer_length'] is not None else 'N/A'
            pen_str = (f"  penalty={prow['primer_penalty']:.2f}"
                       if prow.get('primer_penalty') is not None else '')
            _print(f"  {prow['primer_name']:25s}  "
                  f"Tm={tm_str}  GC={gc_str}  len={len_str}{pen_str}  [{status}]")
        if n_outer_warn > 0:
            _print(f"  Outer primer warnings : {n_outer_warn}, "
                  f"see outer primers CSV for details")
        else:
            _print(f"  Outer primer status   : OK, all primers passed all checks")

    _print("=" * 60)

    if output_csv:
        df.to_csv(output_csv, index=False)
        _print(f"  Oligos CSV saved to   : {output_csv}")
        if not fragment_df.empty:
            frag_csv = output_csv.replace('.csv', '_fragments.csv')
            if not frag_csv.endswith('_fragments.csv'):
                frag_csv = output_csv + '_fragments.csv'
            fragment_df.to_csv(frag_csv, index=False)
            _print(f"  Fragment CSV saved to : {frag_csv}")
        if not outer_df.empty:
            outer_csv = output_csv.replace('.csv', '_outer_primers.csv')
            if not outer_csv.endswith('_outer_primers.csv'):
                outer_csv = output_csv + '_outer_primers.csv'
            outer_df.to_csv(outer_csv, index=False)
            _print(f"  Outer primers saved to: {outer_csv}")
    _print(f"  Runtime               : {_fmt_elapsed(time.time() - _t_start)}")
    _print("=" * 60)

    _print()
    _write_log()
    return df, fragment_df, outer_df


# SECTION 8: EXTERNAL OLIGO EVALUATOR

def evaluate_oligos_csv(
    csv_input,
    main_sequence: str,
    # Quality thresholds -- identical defaults to design_oligos for apples-to-apples comparison
    overlap_tm_range:           Tuple[float, float] = DEFAULT_OVERLAP_TM_RANGE,
    tm_warning_buffer:          float               = DEFAULT_TM_WARNING_BUFFER,
    gc_warning_buffer:          float               = DEFAULT_GC_WARNING_BUFFER,
    overlap_gc_range:           Tuple[float, float] = DEFAULT_OVERLAP_GC_RANGE,
    max_homopolymer_run:        int                 = DEFAULT_MAX_HOMOPOLYMER_RUN,
    max_dinucleotide_repeat:    int                 = DEFAULT_MAX_DINUCLEOTIDE_REPEAT,
    repeat_vicinity_buffer:     int                 = DEFAULT_REPEAT_VICINITY_BUFFER,
    hairpin_dg_threshold:       float               = DEFAULT_HAIRPIN_DG_THRESHOLD,
    hairpin_dg_buffer:          float               = DEFAULT_HAIRPIN_DG_BUFFER,
    homodimer_dg_threshold:     float               = DEFAULT_HOMODIMER_DG_THRESHOLD,
    homodimer_dg_buffer:        float               = DEFAULT_HOMODIMER_DG_BUFFER,
    gc_clamp_length:            int                 = DEFAULT_GC_CLAMP_LENGTH,
    gc_clamp_min_gc:            int                 = DEFAULT_GC_CLAMP_MIN_GC,
    heterodimer_dg_threshold:   float               = DEFAULT_HETERODIMER_DG_THRESHOLD,
    heterodimer_dg_buffer:      float               = DEFAULT_HETERODIMER_DG_BUFFER,
    heterodimer_neighbor_range: int                 = DEFAULT_HETERODIMER_NEIGHBOR_RANGE,
    misprime_length:            int                 = DEFAULT_MISPRIME_LENGTH,
    misprime_tip_length:        int                 = DEFAULT_MISPRIME_TIP_LENGTH,
    misprime_max_mismatches:    int                 = DEFAULT_MISPRIME_MAX_MISMATCHES,
    misprime_kmer_size:         int                 = DEFAULT_MISPRIME_KMER_SIZE,
    # Thermodynamic conditions
    mv_conc:                    float               = DEFAULT_MV_CONC,
    dv_conc:                    float               = DEFAULT_DV_CONC,
    dntp_conc:                  float               = DEFAULT_DNTP_CONC,
    dna_conc:                   float               = DEFAULT_DNA_CONC,
    # Output
    output_csv:                 Optional[str]       = None,
    output_log:                 Optional[str]       = None,
) -> pd.DataFrame:
    """
    Evaluate a set of externally-designed PCA oligos against a target sequence,
    using the same quality checks as design_oligos() for direct benchmarking.

    Parameters
    ----------
    csv_input : str or pd.DataFrame
        CSV file path or DataFrame. Required columns:
            oligo_name     : unique string identifier
            oligo_sequence : 5'->3' ordered oligo sequence (ACGT only)
            overlap_length : int - bases at the 3' end that overlap the next oligo
        Any extra columns will raise a ValueError; only the three required columns
        are accepted. Rows must be in placement order (ascending position along
        main_sequence): the heterodimer check compares each oligo with the oligos
        two rows away.

    main_sequence : str
        Full-length target DNA sequence. Accepts any of:
            - Raw ACGT string
            - A FASTA-formatted string (single or multi-entry; only first entry used)
            - A GenBank-formatted string (whitespace/line-numbers stripped automatically)
            - A file path (.fa, .fasta, .fna, .gb, .gbk, .genbank, etc.) -- the file
              is read and then parsed as FASTA or GenBank as appropriate
        Only unambiguous A, T, G, C bases are allowed. IUPAC ambiguity codes
        (R, Y, S, W, K, M, B, D, H, V, N) will raise a ValueError.

    overlap_tm_range / overlap_gc_range / max_homopolymer_run / ... :
        Quality thresholds. All default to OligoFit's own defaults so the
        comparison is apples-to-apples. Override any parameter to match the
        conditions you care about (e.g. overlap_tm_range=(60, 62)).

    output_csv : str or None
        Save the result DataFrame to this CSV path. Default: None (no file).

    output_log : str or None
        Save a run log to this path. Auto-derived from output_csv when None.
        Set both to None to suppress all disk output. Default: None.

    Returns
    -------
    pd.DataFrame
        One row per oligo with flag columns matching design_oligos output:
        oligo_name, oligo_sequence, overlap_length, overlap_sequence,
        overlap_tm_C, overlap_gc_%, oligo_gc_%, oligo_length,
        direction, global_start, global_end,
        flag_tm_out_of_range, flag_overlap_gc,
        flag_homopolymer, flag_dinucleotide_repeat,
        flag_hairpin, flag_homodimer, flag_gc_clamp,
        flag_heterodimer, flag_misprime,
        flag_body_homopolymer, flag_body_dinucleotide,
        warnings
    """
    import io as _io
    import re as _re
    import sys as _sys
    import time as _time

    _t_start = _time.time()

    # -- Log setup (mirrors design_oligos pattern) --------------------------
    if output_log is None and output_csv is not None:
        _log_stem  = _re.sub(r'\.csv$', '', output_csv, flags=_re.IGNORECASE)
        output_log = _log_stem + '_eval_log.txt'

    _log_buf     = _io.StringIO()
    _real_stdout = _sys.stdout

    def _print(*args, **kwargs) -> None:
        _tmp = _io.StringIO()
        kwargs_copy = {k: v for k, v in kwargs.items() if k != 'file'}
        print(*args, file=_tmp, **kwargs_copy)
        line = _tmp.getvalue()
        _log_buf.write(line)
        try:
            _real_stdout.write(line)
        except UnicodeEncodeError:
            _real_stdout.write(line.encode(_real_stdout.encoding or 'utf-8',
                                           errors='replace').decode(
                                           _real_stdout.encoding or 'utf-8',
                                           errors='replace'))
        _real_stdout.flush()

    _log_written = [False]

    def _write_log() -> None:
        if _log_written[0]:
            return
        _log_written[0] = True
        if output_log:
            try:
                with open(output_log, 'w', encoding='utf-8') as _lf:
                    _lf.write(_log_buf.getvalue())
                _print(f'  Eval log saved to     : {output_log}')
            except OSError as _e:
                _print(f'  [WARNING] Could not write log: {_e}')

    class _LogGuard:
        def __del__(self):
            _write_log()
    _log_guard = _LogGuard()

    def _raise(exc: Exception) -> None:
        _log_buf.write(f'\n[ERROR] {type(exc).__name__}: {exc}\n')
        _write_log()
        raise exc

    # -- Header -------------------------------------------------------------
    _print("=" * 60)
    _print("  evaluate_oligos_csv -- External Oligo Quality Checker")
    _print("=" * 60)
    _print(f"  overlap_tm_range          : {overlap_tm_range}")
    _print(f"  tm_warning_buffer         : {tm_warning_buffer} C")
    _print(f"  gc_warning_buffer         : {gc_warning_buffer} %")
    _print(f"  overlap_gc_range          : {overlap_gc_range}")
    _print(f"  max_homopolymer_run       : {max_homopolymer_run}")
    _print(f"  max_dinucleotide_repeat   : {max_dinucleotide_repeat}")
    _print(f"  repeat_vicinity_buffer    : {repeat_vicinity_buffer}")
    _print(f"  hairpin_dg_threshold      : {hairpin_dg_threshold} kcal/mol")
    _print(f"  hairpin_dg_buffer         : {hairpin_dg_buffer} kcal/mol")
    _print(f"  homodimer_dg_threshold    : {homodimer_dg_threshold} kcal/mol")
    _print(f"  homodimer_dg_buffer       : {homodimer_dg_buffer} kcal/mol")
    _print(f"  gc_clamp_length           : {gc_clamp_length} bp")
    _print(f"  gc_clamp_min_gc           : {gc_clamp_min_gc}")
    _print(f"  heterodimer_dg_threshold  : {heterodimer_dg_threshold} kcal/mol")
    _print(f"  heterodimer_dg_buffer     : {heterodimer_dg_buffer} kcal/mol")
    _print(f"  heterodimer_neighbor_range: {heterodimer_neighbor_range}")
    _print(f"  misprime_length           : {misprime_length}")
    _print(f"  misprime_tip_length       : {misprime_tip_length}")
    _print(f"  misprime_max_mismatches   : {misprime_max_mismatches}")
    _print(f"  misprime_kmer_size        : {misprime_kmer_size}")
    _print(f"  mv_conc={mv_conc} mM  dv_conc={dv_conc} mM  "
           f"dntp_conc={dntp_conc} mM  dna_conc={dna_conc} nM")
    _print("=" * 60)
    _print()

    # -- Step 0: Load CSV / DataFrame ---------------------------------------
    _print("[1/5] Loading and validating input ...")

    if isinstance(csv_input, pd.DataFrame):
        raw_df = csv_input.copy()
    elif isinstance(csv_input, str):
        try:
            raw_df = pd.read_csv(csv_input)
        except FileNotFoundError:
            _raise(FileNotFoundError(f"CSV file not found: {csv_input!r}"))
        except Exception as _e:
            _raise(ValueError(f"Could not read CSV {csv_input!r}: {_e}"))
    else:
        _raise(TypeError(
            f"csv_input must be a file path (str) or pd.DataFrame, "
            f"got {type(csv_input).__name__}"
        ))

    required = {'oligo_name', 'oligo_sequence', 'overlap_length'}
    missing  = required - set(raw_df.columns)
    if missing:
        _raise(ValueError(
            f"CSV is missing required column(s): {sorted(missing)}. "
            f"Required: oligo_name, oligo_sequence, overlap_length"
        ))

    if raw_df.empty:
        _raise(ValueError("CSV contains no rows."))

    extra_cols = [c for c in raw_df.columns if c not in required]
    if extra_cols:
        _raise(ValueError(
            f"CSV contains unexpected column(s): {sorted(extra_cols)}. "
            f"Only 'oligo_name', 'oligo_sequence', and 'overlap_length' are allowed."
        ))

    # -- main_sequence: accept file path, inline FASTA, GenBank, or raw ACGT string
    import os as _os
    _FASTA_EXT  = ('.fa', '.fasta', '.fna', '.ffn', '.frn')
    _GENBANK_EXT = ('.gb', '.gbk', '.genbank', '.gbff')
    _seq_src_label = 'inline'

    _ms_stripped = main_sequence.strip()
    if _ms_stripped and _os.path.isfile(_ms_stripped):
        # It is a real file path -- read its content
        _seq_src_label = _ms_stripped
        try:
            with open(_ms_stripped, 'r', encoding='utf-8', errors='replace') as _sf:
                main_sequence = _sf.read()
        except OSError as _e:
            _raise(OSError(f"Could not read sequence file {_ms_stripped!r}: {_e}"))

    # Parse input sequence (FASTA, GenBank, or raw ACGT)
    main_sequence = clean_input_sequence(main_sequence)

    if not main_sequence:
        _raise(ValueError(
            f"main_sequence is empty after parsing "
            f"({'file: ' + _seq_src_label if _seq_src_label != 'inline' else 'inline string'}). "
            f"Supported formats: raw ACGT string, FASTA, GenBank, or a file path to any of these."
        ))

    _IUPAC_NOTE = (
        "Only unambiguous A, T, G, C bases are supported. "
        "IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N) are not accepted. "
        "Please resolve ambiguous positions before running."
    )
    invalid_chars = set(main_sequence) - set('ACGT')
    if invalid_chars:
        _raise(ValueError(
            f"main_sequence contains invalid character(s): {sorted(invalid_chars)}. "
            + _IUPAC_NOTE
        ))

    _print(f"      {len(raw_df)} oligos loaded  |  main sequence: {len(main_sequence)} bp"
           f"  (source: {_seq_src_label})")

    # Validate each oligo sequence and overlap_length.
    # Terminal oligo (last row) is allowed to have NaN/blank overlap_length, it
    # has no downstream neighbour, so overlap is meaningless. We coerce it to 0.
    # Any non-terminal oligo with a missing or non-integer overlap_length raises.
    last_idx = raw_df.index[-1]
    for idx, row in raw_df.iterrows():
        row_label = idx + 2   # +1 for 0->1-based, +1 for header row visible in the CSV
        seq = str(row['oligo_sequence']).strip().upper()
        if len(seq) == 0:
            _raise(ValueError(
                f"oligo '{row['oligo_name']}' (row {row_label}) has an empty sequence."
            ))
        bad = set(seq) - set('ACGT')
        if bad:
            _raise(ValueError(
                f"oligo '{row['oligo_name']}' (row {row_label}) contains invalid "
                f"character(s): {sorted(bad)}. "
                + _IUPAC_NOTE
            ))
        import math as _math
        raw_val = row['overlap_length']
        is_nan  = (raw_val is None) or (isinstance(raw_val, float) and _math.isnan(raw_val))
        if is_nan:
            if idx == last_idx:
                # Terminal oligo, NaN overlap is expected; treat as 0
                raw_df.at[idx, 'overlap_length'] = 0
                ov = 0
            else:
                _raise(ValueError(
                    f"oligo '{row['oligo_name']}' (row {row_label}) has non-integer "
                    f"overlap_length: {raw_val!r}. Only the terminal oligo may omit this value."
                ))
        else:
            try:
                ov = int(raw_val)
            except (ValueError, TypeError):
                _raise(ValueError(
                    f"oligo '{row['oligo_name']}' (row {row_label}) has non-integer "
                    f"overlap_length: {raw_val!r}"
                ))
        if ov < 0:
            _raise(ValueError(
                f"oligo '{row['oligo_name']}' (row {row_label}) has negative "
                f"overlap_length ({ov})."
            ))
        if ov > len(seq):
            _raise(ValueError(
                f"oligo '{row['oligo_name']}' (row {row_label}) has overlap_length "
                f"({ov}) > oligo length ({len(seq)})."
            ))

    # -- Step 1: Direction resolution ---------------------------------------
    _print("[2/5] Resolving oligo directions against main sequence ...")

    seq_upper = main_sequence
    resolved  = []
    n_fwd, n_rev, n_unres = 0, 0, 0

    for idx, row in raw_df.iterrows():
        oligo_up = str(row['oligo_sequence']).strip().upper()
        name     = str(row['oligo_name'])

        fwd_pos = seq_upper.find(oligo_up)
        if fwd_pos != -1:
            resolved.append({
                'direction':    'forward',
                'global_start': fwd_pos,
                'global_end':   fwd_pos + len(oligo_up),
                'oligo_up':     oligo_up,
                'name':         name,
            })
            n_fwd += 1
            continue

        rc_oligo = reverse_complement(oligo_up)
        rev_pos  = seq_upper.find(rc_oligo)
        if rev_pos != -1:
            resolved.append({
                'direction':    'reverse',
                'global_start': rev_pos,
                'global_end':   rev_pos + len(oligo_up),
                'oligo_up':     oligo_up,
                'name':         name,
            })
            n_rev += 1
            continue

        _warnings.warn(
            f"evaluate_oligos_csv: oligo '{name}' (row {idx}) not found on "
            f"either strand of main_sequence. Direction set to 'unresolvable'. "
            f"Misprime scan will be skipped for this oligo.",
            UserWarning,
            stacklevel=2,
        )
        resolved.append({
            'direction':    'unresolvable',
            'global_start': None,
            'global_end':   None,
            'oligo_up':     oligo_up,
            'name':         name,
        })
        n_unres += 1

    _print(f"      Forward: {n_fwd}  |  Reverse: {n_rev}  |  Unresolvable: {n_unres}")

    # -- Step 2: Per-oligo quality checks -----------------------------------
    _print("[3/5] Running per-oligo quality checks ...")

    per_oligo_records = []

    for i, (raw_idx, row) in enumerate(raw_df.iterrows()):
        res       = resolved[i]
        oligo_up  = res['oligo_up']
        ov_len    = int(row['overlap_length'])
        name      = res['name']
        direction = res['direction']

        warnings_list: List[str] = []
        flags: dict = {
            'flag_tm_out_of_range':     False,
            'flag_overlap_gc':          False,
            'flag_homopolymer':         False,
            'flag_dinucleotide_repeat': False,
            'flag_hairpin':             False,
            'flag_homodimer':           False,
            'flag_gc_clamp':            False,
            'flag_body_homopolymer':    False,
            'flag_body_dinucleotide':   False,
            'flag_heterodimer':         False,
            'flag_misprime':            False,
        }

        # Overlap checks (skip silently when overlap_length == 0)
        overlap_seq = None
        overlap_tm  = None
        overlap_gc  = None
        hp_warn     = []
        dn_warn     = []

        if ov_len > 0:
            # The overlap is always the rightmost N bases of the oligo's top-strand
            # coverage, i.e. fragment[oligo_end - N : oligo_end].
            # Forward oligos are stored as the top strand directly, so the last N bases
            # of the stored sequence is the overlap.
            # Reverse oligos are stored as RC(top_strand_region), so we reconstruct the
            # top strand via RC(oligo_up) and then take its LAST N bases.
            if direction == 'reverse':
                overlap_seq = reverse_complement(oligo_up)[-ov_len:]
            else:
                overlap_seq = oligo_up[-ov_len:]
            overlap_tm  = round(calc_tm(overlap_seq, mv_conc, dv_conc, dntp_conc, dna_conc), 2)
            overlap_gc  = round(gc_percent(overlap_seq), 1)

            min_tm, max_tm = overlap_tm_range
            if not ((min_tm - tm_warning_buffer) <= overlap_tm <= (max_tm + tm_warning_buffer)):
                flags['flag_tm_out_of_range'] = True
                bound   = min_tm if overlap_tm < (min_tm - tm_warning_buffer) else max_tm
                dirword = 'BELOW minimum' if overlap_tm < (min_tm - tm_warning_buffer) else 'ABOVE maximum'
                warnings_list.append(f"Overlap Tm {overlap_tm:.1f}C is {dirword} {bound}C")

            gc_warn = check_overlap_gc(overlap_seq, overlap_gc_range)
            if gc_warn:
                warnings_list.extend(gc_warn)
            # flag rule identical to design_oligos: outside the GC range widened by gc_warning_buffer
            if not ((overlap_gc_range[0] - gc_warning_buffer) <= overlap_gc <= (overlap_gc_range[1] + gc_warning_buffer)):
                flags['flag_overlap_gc'] = True

            hp_warn = check_homopolymer(overlap_seq, max_homopolymer_run)
            if hp_warn:
                flags['flag_homopolymer'] = True
                warnings_list.extend(hp_warn)

            dn_warn = check_dinucleotide_repeat(overlap_seq, max_dinucleotide_repeat)
            if dn_warn:
                flags['flag_dinucleotide_repeat'] = True
                warnings_list.extend(dn_warn)

            # Vicinity check: scanned on main_sequence around the overlap (top-strand coordinates),
            # as in design_oligos, so bases beyond the oligo's end are included
            if direction != 'unresolvable' and ov_len < len(oligo_up):
                overlap_end_g = res['global_end']
                vicinity_warn = check_repeat_vicinity(
                    main_sequence, overlap_end_g - ov_len, overlap_end_g,
                    repeat_vicinity_buffer, max_homopolymer_run, max_dinucleotide_repeat,
                    overlap_hp_results=hp_warn,
                    overlap_dn_results=dn_warn,
                )
                warnings_list.extend(vicinity_warn)

        # Body checks (non-overlap portion): the overlap is the last ov_len bases of a forward oligo
        # and the first ov_len bases of a reverse oligo (stored as the reverse complement)
        if ov_len < len(oligo_up):
            body_seq = oligo_up[ov_len:] if direction == 'reverse' else oligo_up[: len(oligo_up) - ov_len]
            body_hp  = check_homopolymer(body_seq, max_homopolymer_run)
            body_dn  = check_dinucleotide_repeat(body_seq, max_dinucleotide_repeat)
            if body_hp:
                flags['flag_body_homopolymer'] = True
            if body_dn:
                flags['flag_body_dinucleotide'] = True

        # Hairpin
        hp_warn_full, hp_sev = _hairpin_check_sliding(
            oligo_up, hairpin_dg_threshold, mv_conc, dv_conc, dntp_conc, dna_conc,
        )
        # flag and warning both require severity > buffer (same rule as design_oligos)
        if hp_warn_full and hp_sev > hairpin_dg_buffer:
            flags['flag_hairpin'] = True
            warnings_list.extend(hp_warn_full)

        # Homodimer
        hd_warn_full, hd_sev = _homodimer_check_sliding(
            oligo_up, homodimer_dg_threshold, mv_conc, dv_conc, dntp_conc, dna_conc,
        )
        if hd_warn_full and hd_sev > homodimer_dg_buffer:
            flags['flag_homodimer'] = True
            warnings_list.extend(hd_warn_full)

        # GC clamp
        gc_clamp_warn = check_gc_clamp(oligo_up, gc_clamp_min_gc, gc_clamp_length)
        if gc_clamp_warn:
            flags['flag_gc_clamp'] = True
            warnings_list.extend(gc_clamp_warn)

        rec = {
            'oligo_name':      name,
            'oligo_sequence':  oligo_up,
            'overlap_length':  ov_len,
            'overlap_sequence': overlap_seq if overlap_seq else 'N/A',
            'overlap_tm_C':    overlap_tm,
            'overlap_gc_%':    overlap_gc,
            'oligo_gc_%':      round(gc_percent(oligo_up), 1),
            'oligo_length':    len(oligo_up),
            'direction':       direction,
            'global_start':    res['global_start'],
            'global_end':      res['global_end'],
            **flags,
            'warnings':        '; '.join(warnings_list) if warnings_list else 'OK',
            # internal fields removed before output
            '_oligo_up':   oligo_up,
            '_frag_start': res['global_start'],
            '_frag_end':   res['global_end'],
        }

        for col in extra_cols:
            rec[col] = row[col]

        per_oligo_records.append(rec)

    # -- Step 3: Misprime scan ----------------------------------------------
    _print(f"[4/5] Running misprime scan and heterodimer check ...")

    resolvable_indices   = [i for i, r in enumerate(per_oligo_records)
                            if r['_frag_start'] is not None]
    unresolvable_indices = [i for i, r in enumerate(per_oligo_records)
                            if r['_frag_start'] is None]

    if resolvable_indices:
        misprime_records = [
            {
                'oligo_sequence': per_oligo_records[i]['_oligo_up'],
                'fragment':       1,
                'frag_start':     per_oligo_records[i]['_frag_start'],
                'frag_end':       per_oligo_records[i]['_frag_end'],
            }
            for i in resolvable_indices
        ]
        mp_flags, mp_warns = check_misprime_all(
            records            = misprime_records,
            fragments          = [main_sequence],
            misprime_length    = misprime_length,
            tip_length         = misprime_tip_length,
            max_mismatches     = misprime_max_mismatches,
            frag_global_starts = [0],
        )
        for list_pos, orig_idx in enumerate(resolvable_indices):
            rec = per_oligo_records[orig_idx]
            rec['flag_misprime'] = mp_flags[list_pos]
            if mp_warns[list_pos]:
                existing = rec['warnings']
                rec['warnings'] = (mp_warns[list_pos] if existing == 'OK'
                                   else existing + '; ' + mp_warns[list_pos])

    for orig_idx in unresolvable_indices:
        per_oligo_records[orig_idx]['flag_misprime'] = False

    # -- Step 4: Heterodimer check ------------------------------------------
    hetdimer_records = [
        {
            'oligo_sequence': rec['_oligo_up'],
            'fragment':       1,
            'oligo_name':     rec['oligo_name'],
        }
        for rec in per_oligo_records
    ]

    hd_flags, hd_warns = check_heterodimers_all(
        records        = hetdimer_records,
        dg_threshold   = heterodimer_dg_threshold,
        neighbor_range = heterodimer_neighbor_range,
        dg_buffer      = heterodimer_dg_buffer,
        mv_conc        = mv_conc,
        dv_conc        = dv_conc,
        dntp_conc      = dntp_conc,
        dna_conc       = dna_conc,
    )

    for i, rec in enumerate(per_oligo_records):
        rec['flag_heterodimer'] = hd_flags[i]
        if hd_warns[i]:
            existing = rec['warnings']
            rec['warnings'] = (hd_warns[i] if existing == 'OK'
                               else existing + '; ' + hd_warns[i])

    # -- Step 5: Build output DataFrame -------------------------------------
    for rec in per_oligo_records:
        rec.pop('_oligo_up',   None)
        rec.pop('_frag_start', None)
        rec.pop('_frag_end',   None)

    col_order = [
        'oligo_name', 'oligo_sequence', 'overlap_length', 'overlap_sequence',
        'overlap_tm_C', 'overlap_gc_%', 'oligo_gc_%', 'oligo_length',
        'direction', 'global_start', 'global_end',
        'flag_tm_out_of_range', 'flag_overlap_gc',
        'flag_homopolymer', 'flag_dinucleotide_repeat',
        'flag_hairpin', 'flag_homodimer', 'flag_gc_clamp',
        'flag_heterodimer', 'flag_misprime',
        'flag_body_homopolymer', 'flag_body_dinucleotide',
        'warnings',
    ] + extra_cols

    out_df = pd.DataFrame(per_oligo_records)[col_order]

    # -- Summary ------------------------------------------------------------
    n_total    = len(out_df)
    n_warnings = (out_df['warnings'] != 'OK').sum()
    valid_tms  = out_df['overlap_tm_C'].dropna()

    flag_summary = [
        ('flag_tm_out_of_range',     'Tm out of range      '),
        ('flag_overlap_gc',          'Overlap GC out range '),
        ('flag_homopolymer',         'Homopolymer (overlap)'),
        ('flag_dinucleotide_repeat', 'Dinucleotide (overlap)'),
        ('flag_hairpin',             'Potential hairpins   '),
        ('flag_homodimer',           'Homodimer risk       '),
        ('flag_gc_clamp',            'Weak GC clamp        '),
        ('flag_heterodimer',         'Heterodimer risk     '),
        ('flag_misprime',            'Potential misprimes  '),
        ('flag_body_homopolymer',    'Homopolymer (body)   '),
        ('flag_body_dinucleotide',   'Dinucleotide (body)  '),
    ]

    _print("\n" + "=" * 60)
    _print("  evaluate_oligos_csv SUMMARY")
    _print("=" * 60)
    _print(f"  Total oligos          : {n_total}")
    _print(f"  Direction resolved    : {n_fwd + n_rev} "
           f"(fwd={n_fwd}, rev={n_rev}, unresolvable={n_unres})")
    if len(valid_tms):
        _print(f"  Overlap Tm range      : {valid_tms.min():.1f}C -- "
               f"{valid_tms.max():.1f}C  (mean {valid_tms.mean():.1f}C)")
    _print("  --- Quality flag counts ---")
    for col, label in flag_summary:
        _print(f"  {label}: {int(out_df[col].sum())} oligo(s)")
    _print(f"  Total with warnings   : {n_warnings} oligo(s)")
    if n_warnings == 0:
        _print("  Status                : OK -- all oligos passed all checks")
    else:
        _print("  Status                : See 'warnings' column for details")

    # NOTES section -- actionable recommendations based on flag analysis
    if n_total > 0:
        _flag_cols_eval = [col for col, _ in flag_summary]
        _existing_eval = [c for c in _flag_cols_eval if c in out_df.columns]
        _n_flagged_eval = int((out_df[_existing_eval].any(axis=1)).sum())
        _flagged_ratio_eval = _n_flagged_eval / n_total

        _flag_counts_eval = {
            col: int(out_df[col].sum())
            for col in _existing_eval
            if int(out_df[col].sum()) > 0
        }

        _flag_display_eval = {
            'flag_tm_out_of_range':     'Tm out of range',
            'flag_overlap_gc':          'Overlap GC out of range',
            'flag_homopolymer':         'Homopolymer (overlap)',
            'flag_dinucleotide_repeat': 'Dinucleotide repeat (overlap)',
            'flag_hairpin':             'Hairpin',
            'flag_homodimer':           'Homodimer',
            'flag_gc_clamp':            'Weak GC clamp',
            'flag_heterodimer':         'Heterodimer',
            'flag_misprime':            'Misprime',
            'flag_body_homopolymer':    'Homopolymer (body)',
            'flag_body_dinucleotide':   'Dinucleotide repeat (body)',
        }

        if _flagged_ratio_eval >= 0.80 or _flag_counts_eval:
            _print(f"\n  --- NOTES ---")

        if _flagged_ratio_eval >= 0.80:
            _pct_eval = round(_flagged_ratio_eval * 100, 1)
            _print(f"  [!] {_pct_eval}% of oligos ({_n_flagged_eval}/{n_total}) carry "
                   f"at least one flag.")
            _print(f"      Consider re-designing with expanded oligo and overlap length "
                   f"ranges so the")
            _print(f"      algorithm can explore more placement options and find "
                   f"lower-penalty oligo configurations.")

        if _flag_counts_eval:
            _max_cnt_eval = max(_flag_counts_eval.values())
            _top_eval = [col for col, cnt in _flag_counts_eval.items() if cnt == _max_cnt_eval]
            _top_eval_names = [_flag_display_eval.get(f, f) for f in _top_eval]

            if len(_top_eval_names) == 1:
                _print(f"  [i] The most prevalent flag is '{_top_eval_names[0]}' "
                       f"({_max_cnt_eval} oligo(s)).")
            else:
                _joined_eval = "' and '".join(_top_eval_names)
                _print(f"  [i] The most prevalent flags are '{_joined_eval}' "
                       f"({_max_cnt_eval} oligo(s) each).")
            _print(f"      Consider raising its penalty weight and/or expanding the oligo "
                   f"and overlap")
            _print(f"      length ranges so the algorithm has a wider search space to find")
            _print(f"      placements that avoid this issue.")

    _print(f"  Runtime               : {_fmt_elapsed(_time.time() - _t_start)}")
    _print("=" * 60)
    _print()

    if output_csv:
        out_df.to_csv(output_csv, index=False)
        _print(f"  Results saved to      : {output_csv}")

    _write_log()
    return out_df
