import os
import json
import traceback
from fastapi import FastAPI, HTTPException, BackgroundTasks, Header, Depends
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import firebase_admin
from firebase_admin import credentials, firestore
from datetime import datetime, timedelta, timezone

# Import directly since they are now in the same folder
from oligofit import design_oligos, evaluate_dna_sequence, evaluate_oligos_csv

# ---------------------------------------------------------
# Firebase Admin Initialization
# ---------------------------------------------------------
firebase_json_str = os.environ.get("FIREBASE_KEY_JSON")
db = None
if firebase_json_str:
    try:
        cred_dict = json.loads(firebase_json_str)
        cred = credentials.Certificate(cred_dict)
        firebase_admin.initialize_app(cred)
        db = firestore.client()
        print("✅ Firebase Admin initialized successfully.")
    except Exception as e:
        print(f"❌ Failed to initialize Firebase Admin: {e}")
else:
    print("⚠️ FIREBASE_KEY_JSON environment variable not found. Firestore writing will be disabled.")

# ---------------------------------------------------------
# API Key Authentication
# ---------------------------------------------------------
VALID_API_KEYS: set = set()
_api_keys_str = os.environ.get("OLIGOFIT_API_KEYS", "")
if _api_keys_str:
    VALID_API_KEYS = {k.strip() for k in _api_keys_str.split(",") if k.strip()}
    print(f"\U0001f511 Loaded {len(VALID_API_KEYS)} API key(s).")
else:
    print("❌ OLIGOFIT_API_KEYS not set. All API requests will be rejected until keys are configured.")


async def verify_api_key(x_api_key: str = Header(None)):
    """Validate the X-API-Key header against the configured key set."""
    if not VALID_API_KEYS:
        raise HTTPException(
            status_code=503,
            detail="API keys are not configured on the server. Set the OLIGOFIT_API_KEYS environment variable.",
        )
    if not x_api_key or x_api_key not in VALID_API_KEYS:
        raise HTTPException(
            status_code=401,
            detail="Invalid or missing API key. Include your key in the 'X-API-Key' header.",
        )


app = FastAPI(title="OligoFit API")

# Allow CORS for local testing and frontend connection
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Adjust in production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/health")
def health_check():
    """Unauthenticated health check for uptime monitoring."""
    return {"status": "ok", "firestore": db is not None}


class PCADesignRequest(BaseModel):
    orderId: str          # NEW: The Firestore document ID to update
    output_csv: str = "DNA_Sequence"
    sequence: str
    mode: str = "custom"  # "custom" | "short_oligos" | "long_oligos"

    min_frag: int = 400
    max_frag: int = 500
    min_oligo: int = 40
    max_oligo: int = 60
    min_overlap: int = 15
    max_overlap: int = 25

    cut_search: int = 30
    cut_score_win: int = 20
    cut_repeat_penalty: float = 1000.0

    min_tm: float = 60.0
    max_tm: float = 65.0
    min_gc: float = 30.0
    max_gc: float = 65.0
    max_hp: int = 4
    max_dinuc: int = 3
    vicinity: int = 5
    hairpin_dg: float = -9.0
    homodimer_dg: float = -9.0
    clamp_len: int = 3
    clamp_min: int = 1
    heterodimer_dg: float = -9.0
    heterodimer_neighbors: int = 2

    misprime_len: int = 20
    misprime_tip: int = 8
    misprime_mismatch: int = 2
    misprime_kmer: int = 12

    hetdimer_kmer_size: int = 12

    p_tm: float = 20.0
    p_gc: float = 8.0
    p_hp: float = 10.0
    p_dn: float = 6.0
    p_hairpin: float = 10.0
    p_homodimer: float = 6.0
    p_clamp: float = 5.0
    p_misprime: float = 6.0
    p_hetdimer_map: float = 8.0
    p_body_hp: float = 3.0
    p_body_dn: float = 3.0
    p_severity_multiplier: float = 1.0

    tm_warning_buffer: float = 1.0
    gc_warning_buffer: float = 2.0
    hairpin_dg_buffer: float = 1.0
    homodimer_dg_buffer: float = 1.0
    heterodimer_dg_buffer: float = 1.0

    junc_len: int = 25
    junc_min_tm: float = 55.0
    junc_max_tm: float = 68.0
    outer_min_len: int = 15
    outer_max_len: int = 35
    outer_min_tm: float = 55.0
    outer_max_tm: float = 68.0

    mv_conc: float = 50.0
    dv_conc: float = 3.0
    dntp_conc: float = 0.8
    dna_conc: float = 200.0

    assembly_method: str = "gibson"
    goldengate_enzyme: str = "BsaI"

    run_oligos: bool = True
    run_primers: bool = True

class EvaluateRequest(BaseModel):
    sequence: str
    max_homopolymer: int = 4
    max_dinucleotide: int = 3

@app.post("/evaluate")
def evaluate_sequence(req: EvaluateRequest, _: None = Depends(verify_api_key)):
    """Evaluate a DNA sequence for structural quality metrics."""
    try:
        results = evaluate_dna_sequence(
            seq=req.sequence,
            max_homopolymer=req.max_homopolymer,
            max_dinucleotide=req.max_dinucleotide
        )
        return {"status": "success", "data": results}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

class EvaluateCSVRequest(BaseModel):
    csv_content: str
    sequence: str
    output_csv: str = "Evaluated_Oligos"

    min_tm: float = 60.0
    max_tm: float = 65.0
    min_gc: float = 30.0
    max_gc: float = 65.0
    max_hp: int = 4
    max_dinuc: int = 3
    vicinity: int = 5
    hairpin_dg: float = -9.0
    homodimer_dg: float = -9.0
    clamp_len: int = 3
    clamp_min: int = 1
    heterodimer_dg: float = -9.0
    heterodimer_neighbors: int = 2

    misprime_len: int = 20
    misprime_tip: int = 8
    misprime_mismatch: int = 2
    misprime_kmer: int = 12

    tm_warning_buffer: float = 1.0
    gc_warning_buffer: float = 2.0
    hairpin_dg_buffer: float = 1.0
    homodimer_dg_buffer: float = 1.0
    heterodimer_dg_buffer: float = 1.0

    mv_conc: float = 50.0
    dv_conc: float = 3.0
    dntp_conc: float = 0.8
    dna_conc: float = 200.0

@app.post("/evaluate_csv")
def evaluate_csv_endpoint(req: EvaluateCSVRequest, _: None = Depends(verify_api_key)):
    """Evaluate an external CSV of oligos against a target sequence."""
    import io as _io
    import pandas as pd
    import tempfile, os as _os
    
    _log_path = None
    try:
        df = pd.read_csv(_io.StringIO(req.csv_content))
        
        _tf = tempfile.NamedTemporaryFile(
            mode='w', suffix='_log.txt', delete=False, encoding='utf-8'
        )
        _log_path = _tf.name
        _tf.close()

        res_df = evaluate_oligos_csv(
            csv_input=df,
            main_sequence=req.sequence,
            overlap_tm_range=(req.min_tm, req.max_tm),
            overlap_gc_range=(req.min_gc, req.max_gc),
            tm_warning_buffer=req.tm_warning_buffer,
            gc_warning_buffer=req.gc_warning_buffer,
            max_homopolymer_run=req.max_hp,
            max_dinucleotide_repeat=req.max_dinuc,
            repeat_vicinity_buffer=req.vicinity,
            hairpin_dg_threshold=req.hairpin_dg,
            hairpin_dg_buffer=req.hairpin_dg_buffer,
            homodimer_dg_threshold=req.homodimer_dg,
            homodimer_dg_buffer=req.homodimer_dg_buffer,
            gc_clamp_length=req.clamp_len,
            gc_clamp_min_gc=req.clamp_min,
            heterodimer_dg_threshold=req.heterodimer_dg,
            heterodimer_dg_buffer=req.heterodimer_dg_buffer,
            heterodimer_neighbor_range=req.heterodimer_neighbors,
            misprime_length=req.misprime_len,
            misprime_tip_length=req.misprime_tip,
            misprime_max_mismatches=req.misprime_mismatch,
            misprime_kmer_size=req.misprime_kmer,
            mv_conc=req.mv_conc,
            dv_conc=req.dv_conc,
            dntp_conc=req.dntp_conc,
            dna_conc=req.dna_conc,
            output_csv=None,
            output_log=_log_path
        )
        
        results = []
        if res_df is not None and not res_df.empty:
            results.append({"name": f"{req.output_csv}_evaluated.csv", "content": res_df.to_csv(index=False)})
            
        try:
            with open(_log_path, 'r', encoding='utf-8') as _lf:
                results.append({"name": f"{req.output_csv}_log.txt", "content": _lf.read()})
        except Exception:
            pass

        return {"status": "success", "results": results}

    except Exception as e:
        # Read log if exception occurred
        err_results = []
        if _log_path and _os.path.exists(_log_path):
            try:
                with open(_log_path, 'r', encoding='utf-8') as _lf:
                    err_results.append({"name": f"{req.output_csv}_log.txt", "content": _lf.read()})
            except Exception:
                pass
        
        return {"status": "error", "detail": str(e), "results": err_results}
    finally:
        if _log_path and _os.path.exists(_log_path):
            try:
                _os.remove(_log_path)
            except Exception:
                pass

def process_and_save_to_firestore(req: PCADesignRequest):
    """Background worker that runs the 60-min algorithm and uploads results to Firestore."""
    import tempfile, os as _os
    _log_path = None
    try:
        # Create a temp file for the run log BEFORE calling design_oligos so
        # the path can be passed as output_log and written during the run.
        _tf = tempfile.NamedTemporaryFile(
            mode='w', suffix='_log.txt', delete=False, encoding='utf-8'
        )
        _log_path = _tf.name
        _tf.close()

        df, frag_df, primers = design_oligos(
            sequence=req.sequence,
            mode=req.mode,
            fragment_length_range=(req.min_frag, req.max_frag),
            oligo_length_range=(req.min_oligo, req.max_oligo),
            overlap_length_range=(req.min_overlap, req.max_overlap),
            overlap_tm_range=(req.min_tm, req.max_tm),
            overlap_gc_range=(req.min_gc, req.max_gc),
            tm_warning_buffer=req.tm_warning_buffer,
            gc_warning_buffer=req.gc_warning_buffer,
            max_homopolymer_run=req.max_hp,
            max_dinucleotide_repeat=req.max_dinuc,
            repeat_vicinity_buffer=req.vicinity,
            hairpin_dg_threshold=req.hairpin_dg,
            hairpin_dg_buffer=req.hairpin_dg_buffer,
            homodimer_dg_threshold=req.homodimer_dg,
            homodimer_dg_buffer=req.homodimer_dg_buffer,
            gc_clamp_length=req.clamp_len,
            gc_clamp_min_gc=req.clamp_min,
            heterodimer_dg_threshold=req.heterodimer_dg,
            heterodimer_dg_buffer=req.heterodimer_dg_buffer,
            heterodimer_neighbor_range=req.heterodimer_neighbors,
            misprime_length=req.misprime_len,
            misprime_tip_length=req.misprime_tip,
            misprime_max_mismatches=req.misprime_mismatch,
            misprime_kmer_size=req.misprime_kmer,
            hetdimer_kmer_size=req.hetdimer_kmer_size,
            penalty_tm=req.p_tm,
            penalty_gc=req.p_gc,
            penalty_homopolymer=req.p_hp,
            penalty_dinucleotide=req.p_dn,
            penalty_hairpin=req.p_hairpin,
            penalty_homodimer=req.p_homodimer,
            penalty_gc_clamp=req.p_clamp,
            penalty_misprime=req.p_misprime,
            penalty_hetdimer_map=req.p_hetdimer_map,
            penalty_body_homopolymer=req.p_body_hp,
            penalty_body_dinucleotide=req.p_body_dn,
            penalty_severity_multiplier=req.p_severity_multiplier,
            cut_search_window=req.cut_search,
            cut_score_window=req.cut_score_win,
            cut_score_repeat_penalty=req.cut_repeat_penalty,
            junction_overlap_length=req.junc_len,
            junction_tm_range=(req.junc_min_tm, req.junc_max_tm),
            outer_primer_length_range=(req.outer_min_len, req.outer_max_len),
            outer_primer_tm_range=(req.outer_min_tm, req.outer_max_tm),
            assembly_method=req.assembly_method,
            goldengate_enzyme=req.goldengate_enzyme,
            mv_conc=req.mv_conc,
            dv_conc=req.dv_conc,
            dntp_conc=req.dntp_conc,
            dna_conc=req.dna_conc,
            run_oligos=req.run_oligos,
            run_primers=req.run_primers,
            output_csv=None,
            output_log=_log_path,
        )

        results = []
        if df is not None and not df.empty:
            results.append({"name": f"{req.output_csv}_oligos.csv", "content": df.to_csv(index=False)})
        if frag_df is not None and not frag_df.empty:
            results.append({"name": f"{req.output_csv}_fragments.csv", "content": frag_df.to_csv(index=False)})
        if primers is not None and not primers.empty:
            results.append({"name": f"{req.output_csv}_outer_primers.csv", "content": primers.to_csv(index=False)})
        # Always include the run log as the 4th file
        try:
            with open(_log_path, 'r', encoding='utf-8') as _lf:
                results.append({"name": f"{req.output_csv}_log.txt", "content": _lf.read()})
        except Exception as _le:
            print(f"Could not read log file: {_le}")
        finally:
            try:
                _os.remove(_log_path)
            except Exception:
                pass

        # Write results directly to Firestore
        if db and req.orderId:
            db.collection("orders").document(req.orderId).set({
                "status": "Complete",
                "results": results
            }, merge=True)
            print(f"Successfully saved Order {req.orderId} to Firestore.")
        else:
            print("Finished, but db or orderId is missing. Could not save.")

    except Exception as e:
        traceback.print_exc()
        # The log file was written by _raise() before the exception propagated.
        # Read it now so the frontend receives the full parameter dump + [ERROR] line.
        _error_log_content = ""
        try:
            with open(_log_path, 'r', encoding='utf-8') as _lf:
                _error_log_content = _lf.read()
        except Exception:
            pass
        finally:
            try:
                _os.remove(_log_path)
            except Exception:
                pass
        if db and req.orderId:
            db.collection("orders").document(req.orderId).set({
                "status": "Error",
                "errorMessage": str(e),
                "results": [
                    {"name": f"{req.output_csv}_log.txt",
                     "content": _error_log_content}
                ] if _error_log_content else [],
            }, merge=True)


def sweep_old_orders():
    if not db:
        return
    try:
        cutoff = datetime.now(timezone.utc) - timedelta(hours=24)
        # Fetch orders where createdAt is older than 24 hours
        docs = db.collection("orders").where("createdAt", "<", cutoff).stream()
        
        count = 0
        for doc in docs:
            db.collection("orders").document(doc.id).delete()
            count += 1
            
        if count > 0:
            print(f"🧹 Swept and deleted {count} old orders.")
    except Exception as e:
        print(f"⚠️ Failed to cleanup old orders: {e}")


@app.get("/results/{order_id}")
def get_results(order_id: str, _: None = Depends(verify_api_key)):
    """Retrieve the status and results for a submitted design job."""
    if not db:
        raise HTTPException(status_code=503, detail="Firestore is not available.")
    doc = db.collection("orders").document(order_id).get()
    if not doc.exists:
        raise HTTPException(status_code=404, detail=f"Order '{order_id}' not found.")
    data = doc.to_dict()
    return {
        "orderId": order_id,
        "status": data.get("status", "Processing"),
        "results": data.get("results", []),
        "errorMessage": data.get("errorMessage"),
    }


@app.post("/design_oligos")
def design_endpoint(req: PCADesignRequest, background_tasks: BackgroundTasks, _: None = Depends(verify_api_key)):
    # Create the Firestore document upfront so GET /results works immediately.
    # set(merge=True) keeps any fields the frontend has already written.
    if db and req.orderId:
        db.collection("orders").document(req.orderId).set({
            "status": "Processing",
            "sequenceName": req.output_csv,
            "sequenceLength": len(req.sequence),
        }, merge=True)

    # Add the heavy lifting to the background queue
    background_tasks.add_task(process_and_save_to_firestore, req)
    
    # Add the sweeper to keep the database clean
    background_tasks.add_task(sweep_old_orders)
    
    # Return immediately; the job continues in the background.
    return {"status": "Task Started", "orderId": req.orderId}
