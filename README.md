---
title: OligoFit
emoji: 🧬
colorFrom: blue
colorTo: purple
sdk: docker
app_port: 7860
pinned: false
---
# OligoFit

> **Oligo design for Polymerase Cycling Assembly (PCA), solved as an exact lowest-penalty path search (A\* with h = 0, equivalent to Dijkstra's algorithm).**

[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
[![Dependencies: primer3-py · biopython · pandas](https://img.shields.io/badge/deps-primer3--py%20%7C%20biopython%20%7C%20pandas-orange)](requirements.txt)
[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](LICENSE)

---

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Design Modes and Assembly Methods](#design-modes-and-assembly-methods)
- [How It Works](#how-it-works)
- [Quality Checks](#quality-checks)
- [Evaluators](#evaluators)
- [Parameters Reference](#parameters-reference)
- [Output Files](#output-files)
- [Web API](#web-api)
- [Known Limitations](#known-limitations)
- [References](#references)
- [Citation](#citation)
- [License](#license)

---

## Overview

OligoFit designs the oligonucleotides needed to synthesise a DNA sequence by **Polymerase Cycling Assembly (PCA)**. Given a DNA sequence it:

1. cuts the sequence into fragments at positions with low repeat content,
2. tiles every fragment with alternating forward and reverse oligos whose overlaps are optimised for Tm, GC content and sequence quality,
3. reports quality flags and warnings for every oligo, designs the outer PCR primers for every fragment, and prepares the fragment junctions for **Gibson Assembly**, **stitching PCR** or **Golden Gate Assembly**.

OligoFit also contains two evaluation modules: the **Oligo Evaluator**, which scores an oligo set designed by *any* tool with exactly the same checks, and the **DNA Evaluator**, which profiles the repeat content and GC content of a sequence.

Design principles:

- **Exact search per fragment.** The tiling of each fragment is the lowest-penalty path found by an exhaustive shortest-path search, not a greedy or randomised heuristic. (Fragment cut points are chosen sequentially, so the result is optimal per fragment, not for the whole sequence at once.)
- **Deterministic.** The same inputs always give the same outputs.
- **Thermodynamics from primer3.** Melting temperatures and secondary-structure energies come from `primer3-py` (SantaLucia 1998 nearest-neighbour parameters) at 50 mM Na⁺, 3 mM Mg²⁺, 0.8 mM dNTPs and 200 nM oligo concentration; all of these are adjustable. Hairpin and dimer free energies are evaluated at primer3's default of 37 °C.
- **Fragment-local checks.** Each fragment is assembled in its own tube, so misprime and heterodimer checks look only within the fragment being assembled.

---

## Installation

```bash
git clone https://github.com/NoureldenRihan/OligoFit.git
cd OligoFit
pip install -r requirements.txt
```

**Python:** 3.9 or later is required by the dependencies; OligoFit has been tested with Python 3.12 and 3.13.
**Dependencies:** `primer3-py` calculates melting temperatures and secondary-structure energies, `pandas` holds the result tables, and `biopython` provides sequence handling (reverse complements). `requirements.txt` pins primer3-py 2.3.1, pandas 2.2.3 and Biopython 1.88.

### Google Colab / Kaggle

```python
!pip install primer3-py==2.3.1 biopython==1.88 -q
```

Download `oligofit.py` into the working directory (or upload it) and import it as shown below.

---

## Quick Start

```python
from oligofit import design_oligos

sequence = open("my_gene.fasta").read()   # raw sequence, FASTA text, or GenBank text

oligo_df, fragment_df, primer_df = design_oligos(
    sequence=sequence,
    output_csv="my_gene.csv",
)
```

`design_oligos()` returns three DataFrames and prints a progress log and a summary:

| Return value | Contents |
|---|---|
| `oligo_df` | One row per oligo: sequence to order, position, overlap, quality flags, penalty score, warnings |
| `fragment_df` | One row per fragment, with the junction to the next fragment |
| `primer_df` | Two outer PCR primers per fragment |

With `output_csv="my_gene.csv"` four files are written (see [Output Files](#output-files)): `my_gene.csv`, `my_gene_fragments.csv`, `my_gene_outer_primers.csv` and `my_gene_log.txt`. The default name is `mySequence.csv`.

Input may be a raw A/C/G/T string, FASTA (only the first entry is used, with a warning) or multi-line GenBank text; whitespace and line numbers are removed automatically. IUPAC ambiguity codes are not accepted.

Examples of common variations:

```python
# Inspect how the sequence will be fragmented, without designing any oligos
_, frag_only, _ = design_oligos(sequence, run_oligos=False)

# Tile only fragment 3 (1-based) and skip the outer primers
design_oligos(sequence, fragment_indices=[3], run_primers=False)

# Preset for long oligos
design_oligos(sequence, mode="long_oligos")

# Golden Gate with BsaI
design_oligos(sequence, assembly_method="golden_gate", goldengate_enzyme="BsaI")

# Custom thresholds and penalties
design_oligos(sequence, overlap_tm_range=(58.0, 63.0), penalty_tm=30, penalty_gc_clamp=10)
```

---

## Design Modes and Assembly Methods

### Modes

`mode` selects the structural ranges. Only the three ranges below are set by a preset; every other parameter stays user-controllable.

| `mode` | Fragment length | Oligo length | Overlap length |
|---|---|---|---|
| `"custom"` (default) | `fragment_length_range` | `oligo_length_range` | `overlap_length_range` |
| `"short_oligos"` | 400–500 bp | 40–60 bp | 15–25 bp |
| `"long_oligos"` | 3000–4000 bp | 680–700 bp | 30–40 bp |

In `short_oligos` and `long_oligos` mode the preset ranges replace any values you pass for those three parameters; use `"custom"` to set your own.

### Assembly methods

`assembly_method` controls how the fragments are prepared for joining.

| Method | What OligoFit does |
|---|---|
| `"gibson"` (default) | Adjacent fragments share a junction overlap. The overlap is raised to at least 40 bp, and, if `junction_tm_range` is left at its default, its Tm minimum is relaxed to 48 °C (Gibson runs isothermally at about 50 °C). |
| `"stitching_pcr"` | Adjacent fragments share `junction_overlap_length` bp (default 25) for overlap-extension PCR. The overlap Tm is checked and flagged, not optimised. |
| `"golden_gate"` | Fragment overlaps are set to 0 bp. The 4-bp sequence at each fragment boundary is used as the sticky end, and the enzyme cassette (`GCATAC` + recognition site + spacer) is added to the terminal oligos. |

**Golden Gate details**

- Supported enzymes: `BsaI` (default), `BsmBI`, `BbsI`. The input sequence is scanned for internal recognition sites; if the chosen enzyme cuts inside the sequence, the next enzyme without internal sites is used (with a warning). If all three have internal sites, a `ValueError` is raised.
- Overhangs are the **native** 4-bp sequences at the fragment boundaries, not synthetic sequences, so the assembled product keeps the original sequence at every junction.
- A candidate overhang is accepted only if it is not a palindrome, not a homopolymer, does not begin with the enzyme recognition sequence, and differs by at least 2 bases (`GOLDEN_GATE_MIN_OVERHANG_DISTANCE`) from every earlier overhang and from the reverse complement of every earlier overhang. The lowest-scoring valid cut site in the search window is used.
- Because overhangs must be mutually distinct, a single design supports roughly 13–23 fragments (median about 18 in tests on random sequences). If no valid overhang remains, the error message says so; use larger fragments or Gibson Assembly.

---

## How It Works

1. **Input cleaning and validation.** The sequence is cleaned, and the parameters are checked for consistency (for example, overlap length < oligo length, oligo length ≤ fragment length).

2. **Fragmentation.** If the sequence is longer than the maximum fragment length, it is cut into fragments. For each cut, OligoFit takes the midpoint of the fragment length range and searches ±`cut_search_window` bp (default 30) around it. Each candidate cut position is scored over a 20-bp window as *|GC% − 50|* plus 1000 for every homopolymer run or dinucleotide repeat; the lowest score wins. A terminal fragment shorter than the minimum is merged into the previous one. Fragments share a junction overlap (see assembly methods).

3. **Risk maps.** For each fragment, two k-mer indexes (k = 12 by default) are built: one for 3′-end **misprime** risk and one for reverse-complement **heterodimer** risk. They are fast, fragment-local screens used during the search.

4. **Tiling by exact search.** Each fragment is modelled as a directed acyclic graph. A **node** is a candidate oligo placement; an **edge** is a valid overlap junction between successive oligos. In the implementation, placements whose successor starts at the same position on the same strand are merged into one state `(position, parity)`, so the search runs over those states from a single start state at position 0. Expanding a state evaluates every combination of oligo length and overlap length in the allowed ranges; each candidate receives a **penalty** from the quality checks, and the state with the lowest accumulated penalty is expanded next. The heuristic is set to zero (h(n) = 0), so the search is equivalent to Dijkstra's algorithm. It runs until all reachable states are processed, and the cheapest complete path is recovered by backtracking. A runaway-search guard (`max(200 000, fragment length × oligo-length range × overlap-length range × 5)` iterations) stops contradictory constraints from looping forever.

   Penalties are additive:

   ```
   penalty = base_penalty × (1 + penalty_severity_multiplier × severity)
   ```

   where *severity* is how far a check exceeds its tolerated limit (see the table below).

5. **Post-tiling checks.** After the cheapest path is found, a 3′-probe **misprime scan** and a thermodynamic **heterodimer check** are run on the final oligo set (they need the assembled set, so they cannot be part of the search).

6. **Verification.** Each fragment is re-assembled from its oligos and compared with the intended sequence; after all fragments are done, the full assembly is re-assembled from the fragments by the chosen method and compared with the input. A mismatch raises an error.

7. **Outputs.** CSV tables, outer primers and a log are written.

Fragments are tiled in parallel with `concurrent.futures.ProcessPoolExecutor`; a single fragment, or an environment without multiprocessing, runs sequentially with identical results.

---

## Quality Checks

Defaults shown; all are parameters. "Tolerance" is the buffer parameter listed in the reference.

| Check | Default limit | Penalised when (base penalty) | Severity | Flag raised when |
|---|---|---|---|---|
| Overlap Tm | 60–65 °C, tolerance ±1.0 °C | Tm outside 59–66 °C (20) | °C beyond 59 or 66 | same as penalised |
| Overlap GC | 30–65 %, tolerance ±2.0 points | GC outside 30–65 % (8) | percentage points beyond the limit, divided by 10 | GC outside 28–67 % (the warning text appears outside 30–65 %) |
| Homopolymer (overlap) | run ≥ 4 bases | any run in the overlap (10) | number of runs | same as penalised |
| Dinucleotide repeat (overlap) | unit repeated ≥ 3× (≥ 6 bp) | any repeat in the overlap (6) | number of repeats | same as penalised |
| Hairpin | ΔG ≤ −9.0 kcal/mol, tolerance 1.0 | ΔG < −9.0 (10) | kcal/mol below −9.0 | ΔG < −10.0 |
| Homodimer | ΔG ≤ −9.0 kcal/mol, tolerance 1.0 | ΔG < −9.0 (6) | kcal/mol below −9.0 | ΔG < −10.0 |
| GC clamp | ≥ 1 G/C in the last 3 bases | fewer than 1 (5) | none (flat) | same as penalised |
| Misprime | pre-tiling: 12-mer at the 3′ end; post-tiling: 20-bp 3′ probe, 8-bp exact tip, ≤ 2 mismatches in the other 12 bases | pre-tiling only: the 3′ 12-mer occurs elsewhere in the fragment, either strand (6) | number of hits − 1 | post-tiling probe matches another site on either strand, excluding the oligo's own site ±5 bp |
| Heterodimer | pre-tiling: reverse-complement 12-mer; post-tiling: ΔG ≤ −9.0, tolerance 1.0, against the oligo two positions away (±1 skipped) | pre-tiling only: a reverse-complement 12-mer of the oligo occurs elsewhere in the fragment (8) | number of hits − 1 | post-tiling ΔG ≤ −10.0 kcal/mol |
| Body homopolymer / dinucleotide | same limits, on the non-overlap part of the oligo | any run in the body (3 each) | number of runs | `flag_body_homopolymer`, `flag_body_dinucleotide` |
| Repeat proximity | repeats within ±5 bp of the overlap but not inside it | never | none | none (warning text only) |

Penalty = base × (1 + `penalty_severity_multiplier` × severity), so a violation with severity 1 costs twice its base penalty at the default multiplier of 1.0.

Notes:

- Tm and GC thresholds apply to the **overlap**; hairpin, homodimer and GC clamp apply to the **whole oligo**.
- For oligos up to 60 bp, hairpin and homodimer are computed on the whole oligo. For 61–200 bp, 60-bp windows with a 30-bp step plus a 3′-anchored window are used and the worst window is reported. Above 200 bp, only the first and last 60 bp are checked. The heterodimer check uses the 3′-terminal 60 bp of each oligo.
- Tm is not reliable below about 10 bp; OligoFit warns if the minimum overlap length is below 10.

### What "flagged" means in the log

The summary line **"Total with warnings"** counts the oligos whose `warnings` column is not `OK`. It includes repeat-proximity warnings and overlap-GC warnings inside the tolerance band, and it does **not** include the two body-repeat flags. The individual flag columns count exactly the conditions in the "Flag raised when" column.

---

## Evaluators

Both evaluators are public functions of `oligofit.py` and also run behind the web API.

### Oligo Evaluator: `evaluate_oligos_csv`

Scores an oligo set designed by any tool, with the same checks, thresholds and thermodynamic conditions as OligoFit.

```python
import pandas as pd
from oligofit import evaluate_oligos_csv

oligos = pd.DataFrame({
    "oligo_name":     ["O1", "O2", "O3"],
    "oligo_sequence": ["ATG...", "CAT...", "GGA..."],   # placeholders: use full 5'->3' sequences as ordered
    "overlap_length": [20, 20, 0],                      # bases shared with the next oligo
})
result = evaluate_oligos_csv(oligos, main_sequence=open("target.fasta").read(),
                             output_csv="evaluated.csv")
```

- Required columns are exactly `oligo_name`, `oligo_sequence` and `overlap_length` (the last oligo may leave it blank). `csv_input` may be a file path or a DataFrame.
- **Rows must be in placement order** (ascending position along the target): the heterodimer check compares each oligo with the oligos two rows away.
- Each oligo is located on the target by exact match, on either strand (the first match is used). The overlap is the part shared with the next oligo, counted as the rightmost bases of the oligo's coverage in top-strand coordinates.
- Because other tools do not split the sequence into fragments, the misprime and heterodimer checks run against the **whole target**.
- Output: a table with one row per oligo (flags, warnings), a run log (`..._eval_log.txt`) and a summary including "Total with warnings". There is no penalty score, because that is OligoFit's internal search cost.

### DNA Evaluator: `evaluate_dna_sequence`

Profiles a sequence before any design: length, GC content, and the number and identity of homopolymer runs (default ≥ 4 bases) and dinucleotide repeats (default ≥ 3 repeats). It runs automatically at the start of every design and can be used on its own.

```python
from oligofit import evaluate_dna_sequence
evaluate_dna_sequence(seq, max_homopolymer=4, max_dinucleotide=3)
```

---

## Parameters Reference

All parameters of `design_oligos()`. Defaults are read directly from the code.

#### Input and preset

| Parameter | Default | Description |
|---|---|---|
| `sequence` | *required* | Input DNA sequence: raw A/C/G/T, single-entry FASTA, or multi-line GenBank text (whitespace and line numbers are removed; only the first FASTA entry is used). |
| `mode` | `"custom"` | `"custom"`, `"short_oligos"` or `"long_oligos"`. Presets set only the fragment, oligo and overlap length ranges (see Design Modes). |

#### Structural ranges

| Parameter | Default | Description |
|---|---|---|
| `fragment_length_range` | `(400, 500)` | Target fragment length. None = treat entire sequence as one fragment. Cut points are sequence-aware (avoids homopolymers / GC extremes). |
| `oligo_length_range` | `(40, 60)` | Oligo length range. Shorter oligos tried first (lower synthesis cost). |
| `overlap_length_range` | `(15, 25)` | Overlap region shared between consecutive oligos (hybridises during PCA). |

#### Quality thresholds

| Parameter | Default | Description |
|---|---|---|
| `overlap_tm_range` | `(60.0, 65.0)` | Tm range for overlap regions (SantaLucia 1998 nearest-neighbor, salt-corrected). |
| `tm_warning_buffer` | `1.0` | Tolerance in °C added on each side of overlap_tm_range before an overlap Tm is flagged. |
| `overlap_gc_range` | `(30.0, 65.0)` | GC content range for overlap regions. |
| `gc_warning_buffer` | `2.0` | Tolerance in percentage points added on each side of overlap_gc_range before an overlap GC% is flagged. |
| `max_homopolymer_run` | `4` | Flag homopolymer runs of at least this length in the overlap (penalised in the search) and in the oligo body (soft penalty, reported as `flag_body_homopolymer`). |
| `max_dinucleotide_repeat` | `3` | Flag dinucleotide repeats >= this many repetitions in the overlap. |
| `repeat_vicinity_buffer` | `5` | Extra bp scanned on each side of the overlap for repeats just outside it. Produces a proximity warning only (no flag, no penalty). |
| `hairpin_dg_threshold` | `-9.0` | Flag oligos with intramolecular hairpin dG <= this value (kcal/mol). |
| `hairpin_dg_buffer` | `1.0` | Extra kcal/mol beyond hairpin_dg_threshold that is tolerated before a hairpin is flagged. |
| `homodimer_dg_threshold` | `-9.0` | Flag oligos with homodimer dG <= this value (kcal/mol). |
| `homodimer_dg_buffer` | `1.0` | Extra kcal/mol beyond homodimer_dg_threshold that is tolerated before a homodimer is flagged. |
| `gc_clamp_length` | `3` | Number of 3' terminal bases examined for GC clamp check. |
| `gc_clamp_min_gc` | `1` | Minimum G/C count required in the last gc_clamp_length bases. |
| `heterodimer_dg_threshold` | `-9.0` | Flag oligo pairs with heterodimer dG <= this value (kcal/mol). ±1 neighbors skipped (intentional overlap partners). |
| `heterodimer_dg_buffer` | `1.0` | Extra kcal/mol beyond heterodimer_dg_threshold that is tolerated before a heterodimer is flagged. |
| `heterodimer_neighbor_range` | `2` | Neighbor range for heterodimer checks on each side (±1 always skipped). |

#### Misprime and heterodimer detection

| Parameter | Default | Description |
|---|---|---|
| `misprime_length` | `20` | Length of the 3' probe used for post-tiling misprime scan. |
| `misprime_tip_length` | `8` | Exact-match tip length at probe 3' end. Sites must match this tip before the rest of the probe is evaluated. |
| `misprime_max_mismatches` | `2` | Max mismatches in the non-tip probe region for a site to count as a misprime. |
| `misprime_kmer_size` | `12` | k-mer seed length for the pre-tiling per-fragment misprime risk map. |
| `hetdimer_kmer_size` | `12` | k-mer seed length for the pre-tiling heterodimer RC risk map. |

#### Search penalties

| Parameter | Default | Description |
|---|---|---|
| `penalty_tm` | `20` | A* penalty when overlap Tm is outside overlap_tm_range. |
| `penalty_gc` | `8` | A* penalty when overlap GC% is outside overlap_gc_range. |
| `penalty_homopolymer` | `10` | A* penalty for a homopolymer run in the overlap. |
| `penalty_dinucleotide` | `6` | A* penalty for a dinucleotide repeat in the overlap. |
| `penalty_hairpin` | `10` | A* penalty when the full oligo has a stable hairpin. |
| `penalty_homodimer` | `6` | A* penalty when the full oligo has homodimer risk. |
| `penalty_gc_clamp` | `5` | A* penalty for weak 3' GC content. |
| `penalty_misprime` | `6` | A* penalty per off-target k-mer hit in the pre-tiling misprime risk map. |
| `penalty_hetdimer_map` | `8` | A* penalty per unique RC k-mer hit in the pre-tiling hetdimer risk map. |
| `penalty_body_homopolymer` | `3` | Soft A* steering penalty for homopolymer runs in the oligo body (non-overlap). Flagged as flag_body_homopolymer. Set to 0 to disable steering. |
| `penalty_body_dinucleotide` | `3` | Same as penalty_body_homopolymer but for dinucleotide repeats in the body. Flagged as flag_body_dinucleotide. |
| `penalty_severity_multiplier` | `1.0` | Scales each penalty by how far its threshold is exceeded: `penalty = base × (1 + multiplier × severity)`. `0` = flat base penalties. |

#### Fragment cut scoring

| Parameter | Default | Description |
|---|---|---|
| `cut_search_window` | `30` | bp searched on either side of the target cut position for the lowest-penalty fragment boundary. |
| `cut_score_window` | `20` | Sequence window (bp) examined around each candidate fragment cut point. |
| `cut_score_repeat_penalty` | `1000.0` | Penalty per repeat element in the cut_score_window. |

#### Junctions and outer primers

| Parameter | Default | Description |
|---|---|---|
| `junction_overlap_length` | `25` | Shared overlap between consecutive fragments for stitching PCR (last N bp of frag N == first N bp of frag N+1). |
| `junction_tm_range` | `(55.0, 68.0)` | Tm range for junction overlaps (wider than overlap_tm_range). |
| `outer_primer_length_range` | `(15, 35)` | Outer PCR primer length range. For each fragment end, the length whose Tm is closest to the middle of outer_primer_tm_range is chosen. |
| `outer_primer_tm_range` | `(55.0, 68.0)` | Tm range for outer PCR primers. |

#### Assembly method

| Parameter | Default | Description |
|---|---|---|
| `assembly_method` | `"gibson"` | `"gibson"`, `"stitching_pcr"` or `"golden_gate"` (see Design Modes and Assembly Methods). |
| `goldengate_enzyme` | `"BsaI"` | Type IIS enzyme for Golden Gate: "BsaI" \| "BsmBI" \| "BbsI". Auto-switched if the chosen enzyme has internal recognition sites in the input. |

#### Thermodynamic conditions

| Parameter | Default | Description |
|---|---|---|
| `mv_conc` | `50.0` | Monovalent salt concentration in mM (Na+/K+) for Tm calculations. |
| `dv_conc` | `3.0` | Divalent salt concentration in mM (Mg2+). |
| `dntp_conc` | `0.8` | dNTP concentration in mM (chelates Mg2+, reducing effective free Mg2+). |
| `dna_conc` | `200.0` | Total oligo strand concentration in nM. |

#### Workflow, output and execution

| Parameter | Default | Description |
|---|---|---|
| `run_oligos` | `True` | False = skip tiling, misprime, heterodimer, junction, and primer design. Returns a fragment-only table. |
| `run_primers` | `True` | False = skip outer PCR primer design. Oligo and junction tables still produced. |
| `fragment_indices` | `None` | 1-based list of fragments to tile. None = tile all. |
| `n_workers` | `None` | Parallel worker processes for fragment tiling. None = all cores. 1 = sequential. |
| `output_csv` | `"mySequence.csv"` | Output file path (`NAME.csv`); also writes `NAME_fragments.csv` and `NAME_outer_primers.csv`. `None` = do not save. |
| `output_log` | `None` | Log file path. Default: derived from `output_csv` (`NAME_log.txt`). Captures everything printed by the run (not Python warnings). |


---

## Output Files

With `output_csv="NAME.csv"`, `design_oligos()` writes:

| File | Contents |
|---|---|
| `NAME.csv` | All oligo sequences, their placement and length, quality flags, penalty score and warnings |
| `NAME_fragments.csv` | The fragments, with their sequences and the junctions between them |
| `NAME_outer_primers.csv` | The outer PCR primers for each fragment |
| `NAME_log.txt` | The parameters used, a summary of the run and any warnings |

The fragment and primer files are skipped when empty (`run_oligos=False` or `run_primers=False`).

### Oligo table columns

`global_oligo_num`, `oligo_name` (`FragNN_OligoNN_F/R`), `fragment`, `oligo_in_frag`, `direction`, `frag_start`, `frag_end`, `global_start`, `global_end`, `oligo_length`, `oligo_gc_%`, **`oligo_sequence`** (5′→3′ as it should be ordered; already reverse-complemented for reverse oligos), `overlap_sequence`, `overlap_length`, `overlap_tm_C`, `overlap_gc_%`, eleven flag columns (`flag_tm_out_of_range`, `flag_overlap_gc`, `flag_homopolymer`, `flag_dinucleotide_repeat`, `flag_hairpin`, `flag_homodimer`, `flag_gc_clamp`, `flag_heterodimer`, `flag_misprime`, `flag_body_homopolymer`, `flag_body_dinucleotide`), `penalty_score`, `is_last_in_frag`, `warnings` (`OK` or plain-English messages).

For reverse oligos the overlap sequence is the top-strand sequence shared with the next oligo; it lies at the 5′ end of the stored (reverse-complemented) oligo.

### Fragment table columns

`fragment_num`, `fragment_length_bp`, `global_start`, `global_end`, `fragment_sequence`, then the junction to the next fragment: `junction_to_frag`, `junction_start`, `junction_end`, `junction_length`, `junction_sequence`, two method-specific columns (`fwd_gibson_overlap`/`rev_gibson_overlap`, `fwd_stitch_primer`/`rev_stitch_primer`, or `fwd_junction_boundary`/`rev_junction_boundary`), `junction_tm_C`, `junction_gc_%`, junction flags and `warnings`. For Golden Gate the cassettes are included in `fragment_sequence`.

### Outer primer table columns

`primer_name`, `fragment`, `direction`, `primer_sequence`, `primer_length`, `primer_tm_C`, `primer_gc_%`, `primer_penalty`, `flag_tm_out_of_range`, `flag_gc`, `flag_homopolymer_3p`, `flag_hairpin`, `warnings`. For each fragment end the primer length whose Tm is closest to the middle of `outer_primer_tm_range` is chosen. The outer primers use the overlap GC range, the Tm tolerance, the homopolymer limit and the hairpin threshold given for the design.

---

## Web API

`fastapi_backend.py` exposes OligoFit as a FastAPI service (the Docker image serves it on port 7860 and is deployed as a Hugging Face Space). Its extra dependencies are listed in `requirements-server.txt`.

| Endpoint | Method | Purpose |
|---|---|---|
| `/health` | GET | Unauthenticated health check |
| `/design_oligos` | POST | Starts a design job in the background and returns `{"status": "Task Started", "orderId": ...}` |
| `/results/{order_id}` | GET | Status (`Processing`, `Complete`, `Error`) and result files of a job |
| `/evaluate` | POST | DNA Evaluator |
| `/evaluate_csv` | POST | Oligo Evaluator (`csv_content` and `sequence`, plus optional thresholds) |

- **Authentication:** every endpoint except `/health` requires an `X-API-Key` header whose value is in the comma-separated `OLIGOFIT_API_KEYS` environment variable.
- **Results storage:** design results are written to the Firestore document `orders/{orderId}` (service-account JSON in `FIREBASE_KEY_JSON`); orders older than 24 hours are deleted.
- **Defaults:** the request models use the same defaults as `design_oligos()` and `evaluate_oligos_csv()`.
- **Result files from the server** are named `<name>_oligos.csv`, `<name>_fragments.csv`, `<name>_outer_primers.csv` and `<name>_log.txt` (design) and `<name>_evaluated.csv` and `<name>_log.txt` (evaluation).

```bash
docker build -t oligofit .
docker run -p 7860:7860 -e OLIGOFIT_API_KEYS=mykey -e FIREBASE_KEY_JSON='{...}' oligofit
```

Pushes to `main` are synchronised to the Hugging Face Space by `.github/workflows/sync_to_hub.yml`.

---

## Known Limitations

- **Only A, C, G and T** are accepted; IUPAC ambiguity codes must be resolved first.
- **Optimal per fragment, not for the whole sequence.** Fragment cut points are chosen one after another.
- **Long oligos.** Hairpin and homodimer checks cover only the first and last 60 bp of oligos longer than 200 bp, and the heterodimer check uses only the 3′-terminal 60 bp, so internal structures in very long oligos are not screened.
- **Golden Gate** supports roughly 13–23 fragments per design (see above).
- **Oligo Evaluator** assumes rows are in placement order, and locates repeated oligo sequences by their first match.
- **Parallel tiling** needs multiprocessing support; otherwise it falls back to sequential execution.
- **Predicted quality, not experimental success.** All checks are computational predictions. GC-rich or repeat-rich sequences can leave some positions with no clean overlap; the warnings column shows which.
- **Python warnings** are shown on the console but are not written to the run log.

---

## References

1. Stemmer, W.P.C., Crameri, A., Ha, K.D., Brennan, T.M. & Heyneker, H.L. (1995). Single-step assembly of a gene and entire plasmid from large numbers of oligodeoxyribonucleotides. *Gene*, **164**(1), 49–53.
2. SantaLucia, J. Jr. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *PNAS*, **95**(4), 1460–1465.
3. Untergasser, A. et al. (2012). Primer3 – new capabilities and interfaces. *Nucleic Acids Research*, **40**(15), e115.
4. Kwok, S. et al. (1990). Effects of primer-template mismatches on the polymerase chain reaction: human immunodeficiency virus type 1 model studies. *Nucleic Acids Research*, **18**(4), 999–1005.
5. Hart, P.E., Nilsson, N.J. & Raphael, B. (1968). A formal basis for the heuristic determination of minimum cost paths. *IEEE Transactions on Systems Science and Cybernetics*, **4**(2), 100–107.
6. Gibson, D.G. et al. (2009). Enzymatic assembly of DNA molecules up to several hundred kilobases. *Nature Methods*, **6**, 343–345.
7. Engler, C., Kandzia, R. & Marillonnet, S. (2008). A one pot, one step, precision cloning method with high throughput capability. *PLoS ONE*, **3**(11), e3647.
8. Horton, R.M., Hunt, H.D., Ho, S.N., Pullen, J.K. & Pease, L.R. (1989). Engineering hybrid genes without the use of restriction enzymes: gene splicing by overlap extension. *Gene*, **77**(1), 61–68.

---

## Citation

If you use OligoFit, please cite the accompanying manuscript (reference to be added on publication).

---

## License

OligoFit is licensed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0) licence; see [LICENSE](LICENSE). You may use, modify and share it for non-commercial purposes with appropriate credit. Commercial use requires separate written permission from the authors. The third-party packages OligoFit depends on (primer3-py, Biopython, pandas) remain under their own licences.
