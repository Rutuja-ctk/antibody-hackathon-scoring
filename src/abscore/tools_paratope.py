import os
from typing import Optional, List, Tuple
from Bio.PDB import PDBParser, ShrakeRupley


HYDROPHOBIC_RESIDUES = set(["A", "V", "I", "L", "M", "F", "W", "Y", "C"])

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}


def _three_to_one(resname: str) -> Optional[str]:
    return THREE_TO_ONE.get(resname.upper())


def _get_chain_by_id(model, chain_id: str):
    for ch in model.get_chains():
        if ch.id == chain_id:
            return ch
    return None


def _collect_atoms(chain) -> List:
    atoms = []
    for res in chain:
        if res.id[0] != " ":
            continue
        for atom in res:
            atoms.append(atom)
    return atoms


def _paratope_residues(
    heavy_chain,
    light_chain,
    antigen_chains: List,
    distance_cutoff: float = 5.0,
) -> List:
    """
    Define paratope as heavy + light residues that are within distance_cutoff
    of any antigen atom.
    """
    antigen_atoms = []
    for ch in antigen_chains:
        antigen_atoms.extend(_collect_atoms(ch))

    paratope_residues = []
    for chain in [heavy_chain, light_chain]:
        if chain is None:
            continue
        for res in chain:
            if res.id[0] != " ":
                continue
            atoms = [atom for atom in res]
            if not atoms:
                continue
            close = False
            for atom in atoms:
                for ag_atom in antigen_atoms:
                    if atom - ag_atom <= distance_cutoff:
                        close = True
                        break
                if close:
                    break
            if close:
                paratope_residues.append(res)
    return paratope_residues


def compute_paratope_sasa_and_aggscore(
    complex_pdb_path: str,
    heavy_chain_id: str,
    light_chain_id: str,
    antigen_chain_ids: List[str],
) -> Tuple[Optional[float], Optional[float]]:
    """
    Compute:
      - paratope SASA (Å²) for interface residues on heavy + light chains
      - paratope AggScore like: 100 * fraction of hydrophobic paratope residues
    """

    if not os.path.exists(complex_pdb_path):
        return None, None

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", complex_pdb_path)
    model = next(structure.get_models())

    heavy_chain = _get_chain_by_id(model, heavy_chain_id) if heavy_chain_id else None
    light_chain = _get_chain_by_id(model, light_chain_id) if light_chain_id else None

    antigen_chains = []
    for cid in antigen_chain_ids:
        ch = _get_chain_by_id(model, cid)
        if ch is not None:
            antigen_chains.append(ch)

    if heavy_chain is None or light_chain is None or not antigen_chains:
        return None, None

    sr = ShrakeRupley()
    sr.compute(structure, level="A")

    paratope_residues = _paratope_residues(
        heavy_chain, light_chain, antigen_chains, distance_cutoff=5.0
    )
    if not paratope_residues:
        return None, None

    total_sasa = 0.0
    hydro_count = 0
    total_res = 0

    for res in paratope_residues:
        aa = _three_to_one(res.get_resname())
        if aa is None:
            continue
        total_res += 1
        if aa in HYDROPHOBIC_RESIDUES:
            hydro_count += 1
        for atom in res:
            total_sasa += getattr(atom, "sasa", 0.0)

    if total_res == 0:
        return None, None

    frac_hydro = hydro_count / total_res
    aggscore_like = frac_hydro * 100.0

    return total_sasa, aggscore_like
