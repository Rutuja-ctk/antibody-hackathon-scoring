import argparse
import os
import sys
import zipfile
from typing import List, Tuple

import pandas as pd

# ensure src/ is on path
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(CURRENT_DIR)
SRC_DIR = os.path.join(PROJECT_ROOT, "src")
if SRC_DIR not in sys.path:
    sys.path.append(SRC_DIR)

from run_team_folder import process_team_folder  # type: ignore


def extract_master_zip(zip_path: str, extracted_root: str) -> str:
    """
    Extract TEAM_NAME.zip into extracted_root/TEAM_NAME
    and return the path to TEAM_NAME directory.
    """
    zip_path = os.path.abspath(zip_path)
    base = os.path.splitext(os.path.basename(zip_path))[0]  # TEAM_NAME
    dest_dir = os.path.join(extracted_root, base)

    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir, exist_ok=True)
        print(f"Extracting {zip_path} into {dest_dir}")
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(dest_dir)

    inner_dir = os.path.join(dest_dir, base)
    if os.path.isdir(inner_dir):
        return inner_dir
    return dest_dir


def discover_team_roots(root: str, extracted_root: str) -> List[str]:
    """
    Find top level TEAM_NAME roots under 'root'.

    Handles:
      - TEAM_NAME.zip
      - TEAM_NAME directory
    """
    root = os.path.abspath(root)
    os.makedirs(extracted_root, exist_ok=True)

    team_roots: List[str] = []

    for entry in os.listdir(root):
        full = os.path.join(root, entry)

        if os.path.isdir(full):
            team_roots.append(full)
        elif os.path.isfile(full) and entry.lower().endswith(".zip"):
            team_dir = extract_master_zip(full, extracted_root)
            team_roots.append(team_dir)

    return sorted(team_roots)


def find_challenge_dirs(team_root: str) -> List[Tuple[str, str]]:
    """
    Inside TEAM_NAME/, find TEAM_NAME_Challenge1 and TEAM_NAME_Challenge2.

    Returns list of (challenge_name, challenge_dir).
    """
    team_root = os.path.abspath(team_root)
    team_base = os.path.basename(team_root)

    challenge_dirs: List[Tuple[str, str]] = []

    ch1_name = f"{team_base}_Challenge1"
    ch2_name = f"{team_base}_Challenge2"

    ch1_dir = os.path.join(team_root, ch1_name)
    ch2_dir = os.path.join(team_root, ch2_name)

    if os.path.isdir(ch1_dir):
        challenge_dirs.append(("Challenge1", ch1_dir))

    if os.path.isdir(ch2_dir):
        challenge_dirs.append(("Challenge2", ch2_dir))

    return challenge_dirs


def main():
    parser = argparse.ArgumentParser(
        description="Run scoring on many TEAM_NAME master folders or ZIPs and write one combined CSV and Excel."
    )
    parser.add_argument(
        "--root",
        required=True,
        help="Directory containing TEAM_NAME.zip and/or TEAM_NAME folders.",
    )
    parser.add_argument(
        "--dockq_native_ch1",
        default="",
        help="Path to native PDB for DockQ for Challenge 1. Leave empty to disable DockQ.",
    )
    parser.add_argument(
        "--dockq_script",
        default="",
        help="Path to DockQ script (DockQ.py inside DockQ/src/DockQ). Optional.",
    )
    parser.add_argument(
        "--heavy_chain",
        required=True,
        help="Antibody heavy chain ID in the complex PDB (eg A).",
    )
    parser.add_argument(
        "--light_chain",
        required=True,
        help="Antibody light chain ID in the complex PDB (eg B).",
    )
    parser.add_argument(
        "--antigen_chains",
        required=True,
        help="Comma separated antigen chain IDs (eg C or C,D).",
    )
    parser.add_argument(
        "--output",
        default="all_teams_scores",
        help="Base name for combined output (no extension).",
    )

    args = parser.parse_args()

    root = os.path.abspath(args.root)
    extracted_root = os.path.join(root, "_extracted")
    os.makedirs(extracted_root, exist_ok=True)

    team_roots = discover_team_roots(root, extracted_root)
    if not team_roots:
        raise RuntimeError(f"No TEAM_NAME folders or ZIPs found in {root}")

    print("Found master team roots:")
    for tr in team_roots:
        print("  ", tr)

    all_dfs = []

    for team_root in team_roots:
        team_name = os.path.basename(team_root)
        print(f"\n=== Processing team: {team_name} ===")

        challenge_dirs = find_challenge_dirs(team_root)
        if not challenge_dirs:
            print(f"  Warning: no Challenge1 or Challenge2 folders found inside {team_root}")
            continue

        for challenge_name, ch_dir in challenge_dirs:
            print(f"  - Scoring {challenge_name} in {ch_dir}")

            dockq_native = args.dockq_native_ch1 if challenge_name == "Challenge1" else ""

            try:
                df_ch = process_team_folder(
                    team_dir=ch_dir,
                    dockq_native=dockq_native,
                    dockq_script=args.dockq_script or "",
                    heavy_chain=args.heavy_chain,
                    light_chain=args.light_chain,
                    antigen_chains_str=args.antigen_chains,
                )
                df_ch["team_master"] = team_name
                df_ch["challenge"] = challenge_name
                all_dfs.append(df_ch)
            except Exception as e:
                print(f"    Error processing {ch_dir}: {e}")

    if not all_dfs:
        raise RuntimeError("No challenge folders were successfully scored.")

    combined = pd.concat(all_dfs, ignore_index=True)

    out_csv = os.path.join(root, args.output + ".csv")
    out_xlsx = os.path.join(root, args.output + ".xlsx")

    combined.to_csv(out_csv, index=False)
    try:
        combined.to_excel(out_xlsx, index=False)
    except Exception as e:
        print(f"Could not write combined Excel file: {e}")

    print("\nCombined results written to:")
    print("  ", out_csv)
    print("  ", out_xlsx)


if __name__ == "__main__":
    main()
