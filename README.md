
# Antibody Hackathon Scoring Pipeline

This repo provides a small scoring pipeline for antibody hackathon submissions.

For each design it takes:
- One complex structure (antibody + antigen)
- One AlphaFold style PAE JSON
- One FASTA with heavy, light and antigen sequences

and computes:
- Binding and interface quality: PRODIGY, ipSAE, DockQ (Challenge 1 only)
- Developability: NetSolP solubility, CDR SASA proxy
- Novelty: CDRH3 identity vs reference

The outputs are CSV and Excel tables with raw metrics, per metric scores and a final score out of 100.

---

## 1. Submission folder layout

Each team has a master folder named after the team, with one folder per challenge.

```text
TEAM_NAME/
  TEAM_NAME_Challenge1/
    structures/
      design_X_complex.pdb
      design_X_pae.json
    sequences/
      design_X.fasta
    metrics/        (created by the pipeline)
    docs/           (optional)

  TEAM_NAME_Challenge2/
    structures/
      design_Y_complex.pdb
      design_Y_pae.json
    sequences/
      design_Y.fasta
    metrics/
    docs/
````

FASTA format:

```fasta
>Heavy_Chain
HEAVY_CHAIN_SEQUENCE
>Light_Chain
LIGHT_CHAIN_SEQUENCE
>Antigen
ANTIGEN_SEQUENCE
```

The names `design_X` and `design_Y` can be anything. Only extensions matter.

There is a `TEAM_Example` folder in this repo that shows a complete example and can be used to test the pipeline.
There is also a `run_all_teams.py` script and command line to run all team designs in one go (meant for organizers).

---

## 2. External tools needed

Install these tools and make sure they work from the same Python environment.

### Python packages

```bash
pip install "numpy<2.0" pandas biopython openpyxl
```

or similar `conda` commands.

### PRODIGY

Used for binding affinity and contacts.

* Install according to PRODIGY docs.
* Check:

```bash
prodigy -h
```

### ipSAE

Used for ipSAE score and interface pLDDT.

* Install according to ipSAE docs: [https://github.com/DunbrackLab/IPSAE](https://github.com/DunbrackLab/IPSAE)
* Check:

```bash
ipsae -h
```

### NetSolP 1.0

Used for solubility.

* Clone or download NetSolP 1.0 and its ONNX models.
* Check that you can run something like:

```bash
python predict.py --help
```

You must tell the pipeline where NetSolP lives.

Open `src/abscore/tools_developability.py` and set the `NETSOLP_PREDICT` variable so it points to your own `predict.py`:

```python
# Example on Linux/macOS
NETSOLP_PREDICT = "/home/user/tools/NetSolP-1.0/PredictionServer/predict.py"

# Example on Windows
NETSOLP_PREDICT = r"C:\Users\User\NetSolP-1.0\PredictionServer\predict.py"
```

Edit that path to match your setup.

### DockQ (Challenge 1 only)

Used to compare predicted complex to a native reference.

```bash
git clone https://github.com/bjornwallner/DockQ.git
cd DockQ
pip install parallelbar "numpy<2.0"
python src/DockQ/DockQ.py -h
```

You will pass the DockQ script path on the command line.

---

## 3. Environment setup

One example with conda:

```bash
conda create -n abscore python=3.10
conda activate abscore
cd /path/to/this/repo
pip install "numpy<2.0" pandas biopython openpyxl
```

Make sure PRODIGY, ipSAE, NetSolP and DockQ all work inside this environment.

---

## 4. Running the pipeline

There are two scripts in `scripts/`:

* `run_team_folder.py` for a single `TEAM_NAME_ChallengeX` folder.
* `run_all_teams.py` to score many teams at once (for organizers).

You must provide chain IDs of the complex:

* heavy chain ID, for example `A`
* light chain ID, for example `B`
* antigen chain IDs, for example `C` or `C,D`

You also need a reference native complex PDB for DockQ in Challenge 1, for example `native_complex.pdb` in the repo root.

---

### 4.1 Single team - Challenge 1 (with DockQ)

#### Linux / macOS (bash)

```bash
conda activate abscore
cd /path/to/repo

python scripts/run_team_folder.py \
  --team_dir "/path/to/TEAM_NAME/TEAM_NAME_Challenge1" \
  --dockq_native "/path/to/native_complex.pdb" \
  --dockq_script "/path/to/DockQ/src/DockQ/DockQ.py" \
  --heavy_chain A \
  --light_chain B \
  --antigen_chains C
```

#### Windows (Command Prompt)

```bat
conda activate abscore
cd C:\path\to\repo

python scripts\run_team_folder.py ^
  --team_dir "C:\path\to\TEAM_NAME\TEAM_NAME_Challenge1" ^
  --dockq_native "C:\path\to\native_complex.pdb" ^
  --dockq_script "C:\path\to\DockQ\src\DockQ\DockQ.py" ^
  --heavy_chain A ^
  --light_chain B ^
  --antigen_chains C
```

Results are written to:

```text
TEAM_NAME_Challenge1/metrics/hackathon_final_scores.csv
TEAM_NAME_Challenge1/metrics/hackathon_final_scores.xlsx
```

---

### 4.2 Single team - Challenge 2 (no DockQ)

#### Linux / macOS (bash)

```bash
conda activate abscore
cd /path/to/repo

python scripts/run_team_folder.py \
  --team_dir "/path/to/TEAM_NAME/TEAM_NAME_Challenge2" \
  --heavy_chain A \
  --light_chain B \
  --antigen_chains C
```

#### Windows (Command Prompt)

```bat
conda activate abscore
cd C:\path\to\repo

python scripts\run_team_folder.py ^
  --team_dir "C:\path\to\TEAM_NAME\TEAM_NAME_Challenge2" ^
  --heavy_chain A ^
  --light_chain B ^
  --antigen_chains C
```

DockQ is automatically skipped for Challenge 2.

#### Run both challenges back to back for a single team (Windows Command Prompt)

```bat
conda activate abscore
cd C:\path\to\repo

python scripts\run_team_folder.py ^
  --team_dir "C:\path\to\TEAM_NAME\TEAM_NAME_Challenge1" ^
  --dockq_native "C:\path\to\native_complex.pdb" ^
  --dockq_script "C:\path\to\DockQ\src\DockQ\DockQ.py" ^
  --heavy_chain A ^
  --light_chain B ^
  --antigen_chains C ^
&& ^
python scripts\run_team_folder.py ^
  --team_dir "C:\path\to\TEAM_NAME\TEAM_NAME_Challenge2" ^
  --heavy_chain A ^
  --light_chain B ^
  --antigen_chains C
```

---

### 4.3 All teams - both challenges

Put all team master folders or zip files inside a root folder, for example:

```text
ALL_TEAMS/
  TEAM_A/
    TEAM_A_Challenge1/...
    TEAM_A_Challenge2/...
  TEAM_B/
    TEAM_B_Challenge1/...
    TEAM_B_Challenge2/...
  TEAM_C.zip   (contains TEAM_C_Challenge1 and TEAM_C_Challenge2 inside)
```

#### Linux / macOS (bash)

```bash
conda activate abscore
cd /path/to/repo

python scripts/run_all_teams.py \
  --root "/path/to/ALL_TEAMS" \
  --dockq_native_ch1 "/path/to/native_complex.pdb" \
  --dockq_script "/path/to/DockQ/src/DockQ/DockQ.py" \
  --heavy_chain A \
  --light_chain B \
  --antigen_chains C \
  --output "hackathon_all_teams"
```

#### Windows (Command Prompt)

```bat
conda activate abscore
cd C:\path\to\repo

python scripts\run_all_teams.py ^
  --root "C:\path\to\ALL_TEAMS" ^
  --dockq_native_ch1 "C:\path\to\native_complex.pdb" ^
  --dockq_script "C:\path\to\DockQ\src\DockQ\DockQ.py" ^
  --heavy_chain A ^
  --light_chain B ^
  --antigen_chains C ^
  --output "hackathon_all_teams"
```

This writes combined results to:

```text
ALL_TEAMS/hackathon_all_teams.csv
ALL_TEAMS/hackathon_all_teams.xlsx
```

Each row contains team id, challenge id, design id, raw metrics, per metric scores, category scores and final score out of 100.

---

## 5. Notes

* If any tool fails for a design, the script prints a warning and sets that metric to `None` instead of stopping the run.
* DockQ is used only for Challenge 1. For Challenge 2 the weight that would have gone to DockQ is put on ipSAE and related metrics.
* `TEAM_Example` can be used as a template for participants.

```


