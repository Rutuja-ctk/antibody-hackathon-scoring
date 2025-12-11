import os
import re
import subprocess
import sys
from typing import Optional, Tuple

# Default DockQ script path
# On your machine this points to your clone.
# Other users can override this with the DOCKQ_SCRIPT environment variable
# or by passing --dockq_script to the CLI.
DEFAULT_DOCKQ_SCRIPT = os.environ.get(
    "DOCKQ_SCRIPT",
    r"C:\Users\RUTUJA\Documents\DockQ\src\DockQ\DockQ.py",
)


def run_prodigy(
    complex_pdb_path: str,
    prodigy_exe: str = "prodigy",
) -> Tuple[Optional[float], Optional[int]]:
    """
    Run PRODIGY on a complex PDB to get:
      - predicted binding free energy (Delta G, kcal per mol)
      - number of intermolecular contacts

    Returns:
        (delta_g, contacts)
        or (None, None) if PRODIGY is not available or output cannot be parsed.
    """
    complex_pdb_path = os.path.abspath(complex_pdb_path)

    if not os.path.exists(complex_pdb_path):
        raise FileNotFoundError(f"Complex PDB not found for PRODIGY: {complex_pdb_path}")

    try:
        env = os.environ.copy()
        env["PYTHONIOENCODING"] = "utf-8"

        result = subprocess.run(
            [prodigy_exe, complex_pdb_path],
            text=True,
            capture_output=True,
            env=env,
            encoding="utf-8",
            errors="replace",
        )
    except FileNotFoundError:
        print("  Warning: prodigy executable not found. Skipping PRODIGY for this design.")
        return None, None
    except Exception as e:
        print(f"  Warning: PRODIGY error: {e}. Skipping PRODIGY for this design.")
        return None, None

    if result.returncode != 0:
        print(
            f"  Warning: PRODIGY returned non-zero exit code {result.returncode}. "
            f"Will still try to parse output."
        )

    # Combine stdout and stderr for parsing
    output = (result.stdout or "") + "\n" + (result.stderr or "")

    # Debug: print first 20 lines
    print("  Debug: PRODIGY output (first 20 lines):")
    for i, line in enumerate(output.split("\n")[:20], 1):
        print(f"    {i:2d}: {line}")

    # Normalize unicode characters
    output = output.replace("−", "-").replace("–", "-").replace("˚", "").replace("°", "")

    delta_g: Optional[float] = None
    contacts: Optional[int] = None

    # Parse contacts - looking for "No. of intermolecular contacts: 245"
    contact_patterns = [
        r"No\.\s*of\s+intermolecular\s+contacts?\s*:\s*(\d+)",
        r"intermolecular\s+contacts?\s*:\s*(\d+)",
    ]

    for pattern in contact_patterns:
        m = re.search(pattern, output, re.IGNORECASE)
        if m:
            try:
                contacts = int(m.group(1))
                print(f"  ✓ Parsed contacts: {contacts}")
                break
            except ValueError:
                continue

    # Parse Delta G - looking for "Predicted binding affinity (kcal.mol-1):    -21.4"
    dg_patterns = [
        r"Predicted\s+binding\s+affinity\s*\([^)]*\)\s*:\s*([-+]?\d+\.?\d*)",
        r"binding\s+affinity[^:]*:\s*([-+]?\d+\.?\d*)",
        r"Delta\s*G[^:]*:\s*([-+]?\d+\.?\d*)",
        r"ΔG[^:]*:\s*([-+]?\d+\.?\d*)",
    ]

    for pattern in dg_patterns:
        m = re.search(pattern, output, re.IGNORECASE)
        if m:
            try:
                delta_g = float(m.group(1))
                print(f"  ✓ Parsed delta_g: {delta_g}")
                break
            except ValueError:
                continue

    # Fallback: manual search for the binding affinity line
    if delta_g is None:
        print("  Warning: Standard patterns failed. Trying manual extraction...")
        for line in output.split("\n"):
            if "binding affinity" in line.lower() and ":" in line:
                print(f"    Found line: {line}")
                after_colon = line.split(":")[1]
                m2 = re.search(r"([-+]?\d+\.?\d*)", after_colon)
                if m2:
                    try:
                        delta_g = float(m2.group(1))
                        print(f"  ✓ Manually parsed delta_g: {delta_g}")
                        break
                    except ValueError:
                        pass

    if delta_g is None or contacts is None:
        print(f"  ✗ Parsing incomplete: delta_g={delta_g}, contacts={contacts}")

    return delta_g, contacts


def run_ipsae(
    pae_json_path: str,
    complex_pdb_path: str,
    prot1_len: Optional[int] = None,
    prot2_len: Optional[int] = None,
    ipsae_exe: str = "ipsae",
) -> Tuple[Optional[float], Optional[float]]:
    """
    Run ipSAE to get:
      - ipSAE score (0 to 1)
      - mean interface pLDDT (from byres file if available, column 6)

    Returns:
        (ipsae_score, interface_plddt)
        or (None, None) if ipSAE is not available or cannot be parsed.
    """
    pae_json_path = os.path.abspath(pae_json_path)
    complex_pdb_path = os.path.abspath(complex_pdb_path)

    if not os.path.exists(pae_json_path):
        raise FileNotFoundError(f"PAE JSON not found for ipSAE: {pae_json_path}")
    if not os.path.exists(complex_pdb_path):
        raise FileNotFoundError(f"Complex PDB not found for ipSAE: {complex_pdb_path}")

    if prot1_len is None or prot2_len is None:
        print("  Warning: no protein lengths passed to ipSAE. Skipping ipSAE for this design.")
        return None, None

    args = [
        ipsae_exe,
        pae_json_path,
        complex_pdb_path,
        str(prot1_len),
        str(prot2_len),
    ]

    try:
        result = subprocess.run(
            args,
            text=True,
            capture_output=True,
            encoding="utf-8",
            errors="replace",
        )
    except FileNotFoundError:
        print("  Warning: ipsae executable not found. Skipping ipSAE for this design.")
        return None, None
    except Exception as e:
        print(f"  Warning: ipSAE error: {e}. Skipping ipSAE for this design.")
        return None, None

    if result.returncode != 0:
        print(f"  Warning: ipSAE returned non-zero exit code {result.returncode}.")

    pdb_dir = os.path.dirname(complex_pdb_path)
    pdb_basename = os.path.splitext(os.path.basename(complex_pdb_path))[0]

    output_file = os.path.join(pdb_dir, f"{pdb_basename}_{prot1_len}_{prot2_len}.txt")
    byres_file = os.path.join(pdb_dir, f"{pdb_basename}_{prot1_len}_{prot2_len}_byres.txt")

    if not os.path.exists(output_file):
        print(f"  Warning: ipSAE output file not found: {output_file}")
        return None, None

    ipsae_val: Optional[float] = None

    # Parse main output file for ipSAE score
    try:
        with open(output_file, "r", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()

        print(f"  Debug: ipSAE main output has {len(lines)} lines")
        for i, line in enumerate(lines[:5], 1):
            print(f"    {i}: {line.rstrip()}")

        # Expected format: Chn1 Chn2 PAE Dist Type ipSAE ...
        scores = []
        for line in lines:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 6:
                if len(parts) > 4 and parts[4].lower() == "max":
                    try:
                        val = float(parts[5])
                        if 0.0 <= val <= 1.0:
                            scores.append(val)
                            print(
                                f"  Debug: Found ipSAE score {val:.4f} "
                                f"in line: {' '.join(parts[:7])}"
                            )
                    except (ValueError, IndexError):
                        continue

        if scores:
            # Your current version uses min here, keeping that behavior
            ipsae_val = min(scores)
            print(f"  ✓ ipSAE score = {ipsae_val:.4f} (min of {len(scores)} interface scores)")
        else:
            print("  ✗ No valid ipSAE scores found in output file")

    except Exception as e:
        print(f"  Warning: could not parse ipSAE output file: {e}")

    # Parse byres file for interface pLDDT (column 6 = AlignRespLDDT)
    iface_plddt: Optional[float] = None

    if os.path.exists(byres_file):
        try:
            with open(byres_file, "r", encoding="utf-8", errors="replace") as f:
                lines = f.readlines()

            plddt_values = []
            for line in lines[1:]:
                if not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        plddt = float(parts[5])
                        if plddt > 0.0 and plddt <= 100.0:
                            plddt_values.append(plddt)
                    except (ValueError, IndexError):
                        continue

            if plddt_values:
                iface_plddt = sum(plddt_values) / len(plddt_values)
                print(f"  ✓ Interface pLDDT = {iface_plddt:.2f} (from {len(plddt_values)} residues)")
            else:
                iface_plddt = 0.0
                print(
                    "  Info: No non-zero pLDDT values found "
                    "(crystal structure or not available). Setting interface pLDDT = 0.0"
                )

        except Exception as e:
            print(f"  Info: could not parse byres file for pLDDT: {e}")
    else:
        print(f"  Info: byres file not found: {byres_file}")

    return ipsae_val, iface_plddt


def compute_cdr_sasa(
    complex_pdb_path: str,
    heavy_chain: str,
    light_chain: str,
) -> Optional[float]:
    """
    Approximate CDR SASA using PRODIGY contact count as a proxy.

    Idea:
      - PRODIGY already counts intermolecular contacts at the interface.
      - More contacts tends to correlate with more buried paratope area.
      - We map: CDR SASA (Å²) ~ contacts * 25.
    """
    try:
        delta_g, contacts = run_prodigy(complex_pdb_path)
    except Exception as e:
        print(f"  Warning: compute_cdr_sasa could not run PRODIGY again: {e}")
        return None

    if contacts is None:
        print("  Warning: compute_cdr_sasa has no contacts, returning None")
        return None

    cdr_sasa = float(contacts) * 25.0
    print(f"  Proxy CDR SASA from contacts ({contacts}) = {cdr_sasa:.1f} Å^2")
    return cdr_sasa


def run_dockq(
    model_complex_pdb: str,
    native_complex_pdb: str,
    dockq_script: Optional[str] = None,
) -> Optional[float]:
    """
    Run DockQ as a plain Python script:

        python <dockq_script> model.pdb native.pdb

    dockq_script can be passed explicitly, or taken from:
      1) function argument
      2) DEFAULT_DOCKQ_SCRIPT (env DOCKQ_SCRIPT or local default)

    Returns:
        dockq_score in [0, 1], or None if anything fails.
    """
    model_complex_pdb = os.path.abspath(model_complex_pdb)
    native_complex_pdb = os.path.abspath(native_complex_pdb)

    if not os.path.exists(model_complex_pdb):
        raise FileNotFoundError(f"DockQ model complex not found: {model_complex_pdb}")
    if not os.path.exists(native_complex_pdb):
        raise FileNotFoundError(f"DockQ native complex not found: {native_complex_pdb}")

    script = dockq_script or DEFAULT_DOCKQ_SCRIPT
    script = os.path.abspath(script)

    if not os.path.exists(script):
        print(f"  Warning: DockQ script not found at {script}. Skipping DockQ.")
        return None

    cmd = [
        sys.executable,
        script,
        model_complex_pdb,
        native_complex_pdb,
    ]

    try:
        result = subprocess.run(
            cmd,
            text=True,
            capture_output=True,
            encoding="utf-8",
            errors="replace",
        )
    except Exception as e:
        print(f"  Warning: DockQ error: {e}. Skipping DockQ for this design.")
        return None

    output = (result.stdout or "") + "\n" + (result.stderr or "")

    print("  Debug: DockQ output (first 10 lines):")
    for i, line in enumerate(output.splitlines()[:10], 1):
        print(f"    {i:2d}: {line}")

    patterns = [
        r"DockQ\s*[=:]\s*([0-9]*\.[0-9]+)",
        r"^DockQ\s+([0-9]*\.[0-9]+)",
    ]

    score: Optional[float] = None

    for pat in patterns:
        m = re.search(pat, output, flags=re.IGNORECASE | re.MULTILINE)
        if m:
            try:
                v = float(m.group(1))
                if 0.0 <= v <= 1.0:
                    score = v
                    break
            except ValueError:
                continue

    if score is None:
        for line in output.splitlines():
            if "DockQ" not in line:
                continue
            tokens = line.replace("=", " ").replace(":", " ").split()
            for token in tokens:
                try:
                    v = float(token)
                except ValueError:
                    continue
                if 0.0 <= v <= 1.0:
                    score = v
                    break
            if score is not None:
                break

    if score is None:
        print("  Warning: could not parse DockQ score from output.")
    else:
        print(f"  Parsed DockQ score = {score:.3f}")

    return score
