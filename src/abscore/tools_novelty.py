from typing import Tuple, Optional
from Bio import SeqIO
import os

# -------------------------------------------------------------------
# Reference CDRH3 sequences
# -------------------------------------------------------------------

# Challenge 1 reference: Pembrolizumab heavy chain CDRH3 (IMGT style)
# This is the native Keytruda CDRH3 and should NOT be changed.
REF_CDRH3_CH1_PEMBROLIZUMAB = "CARRDYRFDMGFDYW"

# Challenge 2 references:
#   - one human germline-like CDRH3
#   - one Nivolumab CDRH3


REF_CDRH3_CH2_GERMLINE = "CAKYDGIYGELDFW"   
REF_CDRH3_CH2_NIVOLUMAB = "CATNDDYW"  


# -------------------------------------------------------------------
# Helpers to get heavy chain sequence from FASTA
# -------------------------------------------------------------------

def _get_heavy_seq_from_fasta(fasta_path: str) -> Optional[str]:
    """
    Get heavy chain amino acid sequence from FASTA.

    Priority:
      1) Look for a record whose header explicitly mentions "heavy"
         (for example >Heavy_Chain, >HEAVY, >VH_Heavy, etc).
      2) If nothing matches, fall back to the longest sequence in the file.

    This matches your challenge format where teams provide:

        >Heavy_Chain
        HEAVY_SEQUENCE
        >Light_Chain
        LIGHT_SEQUENCE
        >Antigen
        ANTIGEN_SEQUENCE
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        return None

    heavy_rec = None

    # 1) Prefer explicit heavy header
    for rec in records:
        header = (rec.id + " " + rec.description).lower()
        if "heavy" in header:
            heavy_rec = rec
            break

    # 2) Fallback: longest sequence
    if heavy_rec is None:
        heavy_rec = max(records, key=lambda r: len(r.seq))

    return str(heavy_rec.seq).upper().replace("*", "")


# -------------------------------------------------------------------
# Identity computation (sliding window over heavy chain)
# -------------------------------------------------------------------

def _best_window_identity(seq: str, ref: str) -> Optional[float]:
    """
    Slide a window of length len(ref) along seq, compute percent identity
    for each window, and return the maximum identity in percent.

    This is robust even if CDRH3 is shifted a bit but has the same motif.
    """
    L = len(ref)
    if L == 0 or len(seq) < L:
        return None

    best = 0.0
    for i in range(len(seq) - L + 1):
        window = seq[i : i + L]
        matches = sum(1 for a, b in zip(window, ref) if a == b)
        ident = 100.0 * matches / L
        if ident > best:
            best = ident
    return best


def _identity_vs_single_ref(fasta_path: str, ref: str) -> Tuple[Optional[float], Optional[bool]]:
    """
    Core helper: heavy sequence → best window identity vs one reference.
    Returns (identity_percent, qc_flag).
    """
    try:
        heavy = _get_heavy_seq_from_fasta(fasta_path)
    except Exception:
        return None, None

    if not heavy:
        return None, None

    identity = _best_window_identity(heavy, ref)
    if identity is None:
        return None, None

    return identity, True


# -------------------------------------------------------------------
# Public API used by the scoring pipeline
# -------------------------------------------------------------------

def compute_cdr3_identity_from_fasta(fasta_path: str) -> Tuple[Optional[float], Optional[bool]]:
    """
    Compute CDRH3 identity with challenge-aware reference choice.

    Automatic behaviour based on path:
      - If the FASTA path contains "challenge2" (case insensitive),
        we treat this as Challenge 2:
          * compute identity vs human germline CDRH3
          * compute identity vs Nivolumab CDRH3
          * take the MAX of the two identities
            (i.e. "how close are you to any known PD-1 / germline template")

      - Otherwise we treat it as Challenge 1:
          * compute identity vs Pembrolizumab CDRH3 only.

    Returns:
      (identity_percent, qc_flag)

      identity_percent: float in 0 to 100, or None if heavy/ref comparison fails
      qc_flag: True if computation worked, None otherwise
    """
    abs_path = os.path.abspath(fasta_path)
    lower_path = abs_path.lower()

    # Challenge 2: use both germline + Nivolumab
    if "challenge2" in lower_path:
        try:
            heavy = _get_heavy_seq_from_fasta(fasta_path)
        except Exception:
            return None, None

        if not heavy:
            return None, None

        id_germ = _best_window_identity(heavy, REF_CDRH3_CH2_GERMLINE)
        id_nivo = _best_window_identity(heavy, REF_CDRH3_CH2_NIVOLUMAB)

        vals = [v for v in (id_germ, id_nivo) if v is not None]
        if not vals:
            return None, None

        identity = max(vals)  # "closest to any of the two references"
        return identity, True

    # Challenge 1: Pembrolizumab only
    return _identity_vs_single_ref(fasta_path, REF_CDRH3_CH1_PEMBROLIZUMAB)


# Backwards-compatible alias used everywhere else
def compute_cdr3_identity(fasta_path: str) -> Tuple[Optional[float], Optional[bool]]:
    """
    Thin wrapper so older code that calls compute_cdr3_identity still works.
    """
    return compute_cdr3_identity_from_fasta(fasta_path)
