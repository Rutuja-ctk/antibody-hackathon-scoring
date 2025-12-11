import argparse
import os
import sys
from typing import List, Tuple

import pandas as pd

# add src/ to path
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(CURRENT_DIR)
SRC_DIR = os.path.join(PROJECT_ROOT, "src")
if SRC_DIR not in sys.path:
    sys.path.append(SRC_DIR)

from abscore.models import DesignInput, RawMetrics, row_from_raw_and_scores  # type: ignore
from abscore.scoring import score_design  # type: ignore
from abscore.tools_binding import (  # type: ignore
    run_prodigy,
    run_ipsae,
    run_dockq,
    compute_cdr_sasa,
)
from abscore.tools_developability import compute_developability_metrics  # type: ignore
from abscore.tools_novelty import compute_cdr3_identity  # type: ignore


def find_design_ids(team_dir: str) -> List[str]:
    """Find all design IDs by scanning structures/ for *_complex.pdb."""
    struct_dir = os.path.join(team_dir, "structures")
    if not os.path.isdir(struct_dir):
        raise FileNotFoundError(f"No structures/ directory found in {team_dir}")

    design_ids: List[str] = []
    for fname in os.listdir(struct_dir):
        if fname.endswith("_complex.pdb"):
            base = fname[: -len("_complex.pdb")]
            design_ids.append(base)

    return sorted(set(design_ids))


def infer_lengths_from_fasta(fasta_path: str) -> Tuple[int, int]:
    """
    Infer antibody length (VH+VL) and antigen length from FASTA.

    Expects headers in the FASTA like:
      >Heavy_Chain
      >Light_Chain
      >Antigen

    Returns:
      (antibody_len, antigen_len)
    """
    antibody_len = 0
    antigen_len = 0

    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")

    label = None
    seq_chunks = {"heavy": [], "light": [], "antigen": []}

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].strip().lower()
                if "heavy" in header:
                    label = "heavy"
                elif "light" in header:
                    label = "light"
                elif "antigen" in header or "ag" in header:
                    label = "antigen"
                else:
                    label = None
            else:
                if label is not None:
                    seq_chunks[label].append(line)

    heavy_seq = "".join(seq_chunks["heavy"])
    light_seq = "".join(seq_chunks["light"])
    antigen_seq = "".join(seq_chunks["antigen"])


    antibody_len = len(heavy_seq) + len(light_seq)
    antigen_len = len(antigen_seq)

    if antibody_len == 0 or antigen_len == 0:
        print(
            f"  Warning: could not infer antibody/antigen lengths from {os.path.basename(fasta_path)}."
        )

    return antibody_len, antigen_len


def process_team_folder(
    team_dir: str,
    dockq_native: str = "",
    dockq_script: str = "",
    heavy_chain: str = "",
    light_chain: str = "",
    antigen_chains_str: str = "",
) -> pd.DataFrame:
    """
    Process one TEAM_NAME_ChallengeX folder and return a DataFrame of scores.

    This is used both by this script and by scripts/run_all_teams.py
    """
    team_dir = os.path.abspath(team_dir)
    team_name = os.path.basename(team_dir)

    print(f"Processing team folder: {team_dir}")

    metrics_dir = os.path.join(team_dir, "metrics")
    os.makedirs(metrics_dir, exist_ok=True)

    design_ids = find_design_ids(team_dir)
    if not design_ids:
        raise RuntimeError(f"No *_complex.pdb designs found in {team_dir}/structures")

    rows = []

    for design_id in design_ids:
        print(f"\n  === Design: {design_id} ===")

        complex_pdb = os.path.join(team_dir, "structures", f"{design_id}_complex.pdb")
        pae_json = os.path.join(team_dir, "structures", f"{design_id}_pae.json")
        fasta_path = os.path.join(team_dir, "sequences", f"{design_id}.fasta")

        missing = [p for p in [complex_pdb, pae_json, fasta_path] if not os.path.exists(p)]
        if missing:
            print("  Missing core files, skipping this design:")
            for m in missing:
                print(f"    - {m}")
            continue

        raw = RawMetrics()

        # 1. PRODIGY: DeltaG and contacts
        try:
            delta_g, contacts = run_prodigy(complex_pdb)
            raw.delta_g = delta_g
            raw.contacts = contacts
            if delta_g is not None:
                print(f"  PRODIGY ΔG = {delta_g:.2f}")
            if contacts is not None:
                print(f"  PRODIGY contacts = {contacts}")
        except Exception as e:
            print(f"  Warning: PRODIGY failed: {e}")

        # 2. ipSAE: ipSAE score and interface pLDDT
        prot1_len, prot2_len = infer_lengths_from_fasta(fasta_path)
        if prot1_len > 0 and prot2_len > 0:
            try:
                ipsae, iface_plddt = run_ipsae(
                    pae_json_path=pae_json,
                    complex_pdb_path=complex_pdb,
                    prot1_len=prot1_len,
                    prot2_len=prot2_len,
                )
                raw.ipsae = ipsae
                raw.iface_plddt = iface_plddt
            except Exception as e:
                print(f"  Warning: ipSAE failed: {e}")
        else:
            print("  Warning: skipping ipSAE because lengths could not be inferred from FASTA.")

        # 3. DockQ: only if dockq_native is provided
        if dockq_native:
            try:
                dockq_val = run_dockq(
                    model_complex_pdb=complex_pdb,
                    native_complex_pdb=dockq_native,
                    dockq_script=dockq_script or None,
                )
                raw.dockq = dockq_val
                if dockq_val is not None:
                    print(f"  DockQ = {dockq_val:.3f}")
            except Exception as e:
                print(f"  Warning: DockQ failed: {e}")
        else:
            raw.dockq = None

        # 4. Developability (NetSolP)
        try:
            netsolp_score, _ = compute_developability_metrics(fasta_path)
            raw.netsolp = netsolp_score
        except Exception as e:
            print(f"  Warning: NetSolP failed: {e}")

        # 5. CDR SASA
        try:
            cdr_sasa_val = compute_cdr_sasa(complex_pdb, heavy_chain, light_chain)
            raw.cdr_sasa = cdr_sasa_val
            if cdr_sasa_val is not None:
                print(f"  CDR SASA = {cdr_sasa_val:.1f} Å^2")
        except Exception as e:
            print(f"  Warning: CDR SASA computation failed: {e}")

        # 6. Novelty (CDR3 identity vs reference)
        try:
            cdr3_id, anarci_pass = compute_cdr3_identity(fasta_path)
            raw.cdr3_identity = cdr3_id
            raw.anarci_pass = anarci_pass
            if cdr3_id is not None:
                print(f"  CDR3 identity = {cdr3_id:.1f} percent")
        except Exception as e:
            print(f"  Warning: CDR3 identity computation failed: {e}")

        # 7. Score everything
        design = DesignInput(team_name=team_name, design_id=design_id)
        scores = score_design(raw)
        row = row_from_raw_and_scores(design, raw, scores)
        rows.append(row)

        print(f"  Viable: {scores.is_viable}, final score: {scores.final_score_100:.1f}")

    if not rows:
        raise RuntimeError("No designs were successfully processed")

    df = pd.DataFrame(rows)

    # write per team outputs
    out_csv = os.path.join(metrics_dir, "hackathon_final_scores.csv")
    out_xlsx = os.path.join(metrics_dir, "hackathon_final_scores.xlsx")
    df.to_csv(out_csv, index=False)
    try:
        df.to_excel(out_xlsx, index=False)
    except Exception as e:
        print(f"Could not write Excel in {metrics_dir}: {e}")

    print("\nWrote scores to:")
    print(" ", out_csv)
    print(" ", out_xlsx)

    return df


def main():
    parser = argparse.ArgumentParser(
        description="Run scoring on a single TEAM_NAME_ChallengeX folder."
    )
    parser.add_argument(
        "--team_dir",
        required=True,
        help="Path to TEAM_NAME_ChallengeX folder (containing structures/ and sequences/).",
    )
    parser.add_argument(
        "--dockq_native",
        default="",
        help="Path to native complex PDB for DockQ. Leave empty to disable DockQ.",
    )
    parser.add_argument(
        "--dockq_script",
        default="",
        help="Path to DockQ script (DockQ.py inside DockQ/src/DockQ). Optional.",
    )
    parser.add_argument(
        "--heavy_chain",
        required=True,
        help="Heavy chain ID in complex (eg A).",
    )
    parser.add_argument(
        "--light_chain",
        required=True,
        help="Light chain ID in complex (eg B).",
    )
    parser.add_argument(
        "--antigen_chains",
        required=True,
        help="Comma separated antigen chain IDs (eg C or C,D).",
    )

    args = parser.parse_args()

    process_team_folder(
        team_dir=args.team_dir,
        dockq_native=args.dockq_native,
        dockq_script=args.dockq_script,
        heavy_chain=args.heavy_chain,
        light_chain=args.light_chain,
        antigen_chains_str=args.antigen_chains,
    )


if __name__ == "__main__":
    main()
