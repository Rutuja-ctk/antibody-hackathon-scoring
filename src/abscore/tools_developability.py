import os
import csv
import subprocess
from typing import Optional, Tuple


def run_netsolp(
    fasta_path: str,
    netsolp_root: str = r"C:\Users\RUTUJA\Documents\NetSolP-1.0",
    python_exe: str = "python",
    model_type: str = "ESM12",
    prediction_type: str = "S",
) -> Optional[float]:
    """
    Run NetSolP on a FASTA file and return a single solubility score
    between 0 and 1 for the Fv.

    This calls:
        cd NetSolP-1.0/PredictionServer
        python predict.py --FASTA_PATH <fasta> --OUTPUT_PATH <csv> \
                          --MODEL_TYPE ESM12 --PREDICTION_TYPE S
    """
    fasta_path = os.path.abspath(fasta_path)
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA not found for NetSolP: {fasta_path}")

    pred_server_dir = os.path.join(netsolp_root, "PredictionServer")
    predict_script = os.path.join(pred_server_dir, "predict.py")

    if not os.path.exists(predict_script):
        print("  Warning: NetSolP predict.py not found at:")
        print(f"           {predict_script}")
        print("           Please update netsolp_root in tools_developability.py")
        return None

    base = os.path.splitext(os.path.basename(fasta_path))[0]
    output_csv = os.path.join(os.path.dirname(fasta_path), f"{base}_netsolp.csv")

    cmd = [
        python_exe,
        predict_script,
        "--FASTA_PATH",
        fasta_path,
        "--OUTPUT_PATH",
        output_csv,
        "--MODEL_TYPE",
        model_type,
        "--PREDICTION_TYPE",
        prediction_type,
    ]

    try:
        result = subprocess.run(
            cmd,
            cwd=pred_server_dir,
            text=True,
            capture_output=True,
            encoding="utf-8",
            errors="replace",
        )
    except Exception as e:
        print(f"  Warning: NetSolP execution error: {e}")
        return None

    if result.returncode != 0:
        print("  Warning: NetSolP returned non zero exit code.")
        print("  stdout:")
        print(result.stdout)
        print("  stderr:")
        print(result.stderr)
        return None

    if not os.path.exists(output_csv):
        print(f"  Warning: NetSolP did not create output file {output_csv}")
        return None

    try:
        with open(output_csv, "r", newline="", encoding="utf-8") as f:
            reader = csv.reader(f)
            header = next(reader, None)
            row = next(reader, None)
            if row is None:
                print("  Warning: NetSolP CSV had no data rows.")
                return None

            best_val = None
            for cell in row:
                try:
                    val = float(cell)
                except ValueError:
                    continue
                if 0.0 <= val <= 1.0:
                    if best_val is None or val > best_val:
                        best_val = val

            if best_val is None:
                print("  Warning: could not find numeric probability in NetSolP output.")
                print("           Row was:", row)
                return None

            print(f"  NetSolP solubility score = {best_val:.3f}")
            return best_val

    except Exception as e:
        print(f"  Warning: error reading NetSolP CSV {output_csv}: {e}")
        return None


def compute_developability_metrics(
    fasta_path: str,
    netsolp_root: str = r"C:\Users\RUTUJA\Documents\NetSolP-1.0",
) -> Tuple[Optional[float], dict]:
    """
    Compute developability related metrics from sequence.

    Returns:
      (netsolp_score, extra_info_dict)
    """
    netsolp_score = run_netsolp(fasta_path, netsolp_root=netsolp_root)
    return netsolp_score, {}
