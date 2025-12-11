from typing import Tuple, Optional
from Bio import SeqIO
import os
import re

# Pembrolizumab heavy chain CDRH3 (IMGT like, includes leading C and trailing W)
REF_CDRH3 = "CARRDYRFDMGFDYW"


def _get_heavy_seq_from_fasta(fasta_path: str) -> Optional[str]:
    """
    Heuristic:
      - Parse all sequences from FASTA
      - Take the longest sequence as heavy chain

    This works for:
      - Full IgG sequences (heavy is longest)
      - VH only constructs
    """
    if not os.path.exists(fasta_path):
        print(f"  Warning: FASTA not found for novelty: {fasta_path}")
        return None

    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        print(f"  Warning: no sequences found in FASTA for novelty: {fasta_path}")
        return None

    rec = max(records, key=lambda r: len(r.seq))
    return str(rec.seq).upper().replace("*", "")


def _extract_cdrh3(seq: str) -> Optional[str]:
    """
    Heuristic CDRH3 extraction using an IMGT style motif.

    We look for the pattern:

        C [5 to 25 any residues] W G Q G

    which matches typical CDRH3 segments like:

        ...VYYCARRDYRFDMGFDYWGQGTTV...

    If that fails, we fall back to:
        - last W in the sequence
        - last C before that W
      and return C...W if length is at least 5.
    """
    # Primary motif based on IMGT like pattern
    m = re.search(r"(C[A-Z]{5,25}W)GQG", seq)
    if m:
        return m.group(1)

    # Fallback: last C before last W
    w_pos = seq.rfind("W")
    if w_pos == -1:
        return None
    prev = seq[:w_pos]
    c_pos = prev.rfind("C")
    if c_pos == -1:
        return None
    if w_pos - c_pos < 5:
        return None

    return seq[c_pos : w_pos + 1]


def compute_cdr3_identity(fasta_path: str) -> Tuple[Optional[float], Optional[bool]]:
    """
    Compute CDRH3 identity vs Pembrolizumab from FASTA only.

    Steps:
      1. Parse FASTA and pick heavy chain as the longest sequence.
      2. Extract CDRH3 from heavy chain using IMGT like motif.
      3. Align CDRH3 to REF_CDRH3 and compute percent identity
         on the overlapping length.

    Returns:
      (identity_percent, qc_flag)

      identity_percent:
        - 0 to 100
        - 0.0 if CDRH3 could not be extracted (so it falls into "novel" band)
      qc_flag:
        - True  if extraction looked OK
        - False if we could not find a clean CDRH3 motif

    This matches your manual:
      - Uses IMGT style CDRH3 motif (C ... W)
      - Uses Pembrolizumab CDRH3 as reference
    """
    try:
        heavy = _get_heavy_seq_from_fasta(fasta_path)
    except Exception as e:
        print(f"  Warning: novelty FASTA parse error: {e}")
        return 0.0, False

    if not heavy:
        print("  Warning: no heavy like sequence found, setting CDR3 identity = 0")
        return 0.0, False

    cdr3 = _extract_cdrh3(heavy)
    if not cdr3:
        print("  Warning: could not extract CDRH3 motif, setting identity = 0")
        return 0.0, False

    ref = REF_CDRH3
    L = min(len(ref), len(cdr3))
    if L == 0:
        print("  Warning: CDRH3 or reference length is zero, setting identity = 0")
        return 0.0, False

    matches = sum(1 for i in range(L) if ref[i] == cdr3[i])
    identity = 100.0 * matches / L

    print(
        f"  CDRH3 identity vs Pembrolizumab: {identity:.1f} percent "
        f"(len(ref)={len(ref)}, len(cdr3)={len(cdr3)})"
    )

    return identity, True
