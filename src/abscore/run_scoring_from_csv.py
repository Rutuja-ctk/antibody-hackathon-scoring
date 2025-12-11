import argparse
import os
import sys
import pandas as pd

# Make sure Python can find src\abscore
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(CURRENT_DIR)
SRC_DIR = os.path.join(PROJECT_ROOT, "src")
if SRC_DIR not in sys.path:
    sys.path.append(SRC_DIR)

from abscore.models import DesignInput, RawMetrics
from abscore.scoring import score_design, build_result_row


def main():
    parser = argparse.ArgumentParser(
        description="Score antibody designs from a CSV of raw metrics."
    )
    parser.add_argument("--input_csv", required=True, help="Input CSV file")
    parser.add_argument("--output_csv", required=True, help="Output CSV file")
    args = parser.parse_args()

    df_in = pd.read_csv(args.input_csv)
    rows = []

    for _, row in df_in.iterrows():
        design = DesignInput(
            team_name=str(row["team_name"]),
            design_id=str(row["design_id"]),
        )

        raw = RawMetrics(
            delta_g=row.get("delta_g"),
            contacts=row.get("contacts"),
            dockq=row.get("dockq"),
            iface_plddt=row.get("iface_plddt"),
            ipsae=row.get("ipsae"),
            camsol=row.get("camsol"),
            aggscore=row.get("aggscore"),
            cdr3_identity=row.get("cdr3_identity"),
            anarci_pass=bool(row.get("anarci_pass")) if "anarci_pass" in row else None,
        )

        scores = score_design(raw)
        out_row = build_result_row(design, raw, scores)
        rows.append(out_row)

    df_out = pd.DataFrame(rows)
    df_out.to_csv(args.output_csv, index=False)
    print(f"Wrote {len(df_out)} scored designs to {args.output_csv}")


if __name__ == "__main__":
    main()
