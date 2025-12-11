import os
import subprocess
from typing import Optional


def run_dockq(
    native_pdb_path: str,
    model_pdb_path: str,
    dockq_exe: str = "DockQ.py",
) -> Optional[float]:
    """
    Run DockQ on a model complex against a native complex.

    Returns:
        DockQ score as float, or None if DockQ is not available or fails.
    """
    native_pdb_path = os.path.abspath(native_pdb_path)
    model_pdb_path = os.path.abspath(model_pdb_path)

    if not os.path.exists(native_pdb_path):
        raise FileNotFoundError(f"DockQ native PDB not found: {native_pdb_path}")
    if not os.path.exists(model_pdb_path):
        raise FileNotFoundError(f"DockQ model PDB not found: {model_pdb_path}")

    try:
        result = subprocess.run(
            [dockq_exe, native_pdb_path, model_pdb_path],
            check=True,
            text=True,
            capture_output=True,
        )
    except FileNotFoundError:
        print("  Warning: DockQ executable not found. Skipping DockQ for this design.")
        return None
    except subprocess.CalledProcessError as e:
        print(f"  Warning: DockQ failed: {e}")
        return None

    dockq_score = None
    for line in result.stdout.splitlines():
        if "DockQ" in line:
            parts = line.replace(":", " ").split()
            for token in parts:
                try:
                    val = float(token)
                    dockq_score = val
                    break
                except ValueError:
                    continue
        if dockq_score is not None:
            break

    if dockq_score is None:
        print("  Warning: could not parse DockQ score from output, setting DockQ=None")
        return None

    return dockq_score
