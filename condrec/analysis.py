from typing import Dict, List, Tuple


def _require_rdkit():
    try:
        from rdkit import Chem  # noqa: F401
        from rdkit.Chem import Descriptors  # noqa: F401
        from rdkit.Chem import rdMolDescriptors  # noqa: F401
    except Exception:
        raise ImportError(
            "RDKit is required. Install via conda: 'conda install -c conda-forge rdkit' "
            "or pip (platform-dependent): 'pip install rdkit-pypi'"
        )


def _split_reaction_smiles(rxn: str) -> Tuple[List[str], List[str], List[str]]:
    parts = rxn.strip().split(">")
    if len(parts) not in (2, 3):
        # Normalize to 3 parts: reactants>agents>products
        # If only reactants>products, insert empty agents
        if len(parts) == 2:
            parts = [parts[0], "", parts[1]]
        else:
            raise ValueError("Reaction SMILES must have 2 or 3 parts separated by '>'")
    reactants, agents, products = parts
    to_list = lambda s: [p for p in s.split(".") if p] if s else []
    return to_list(reactants), to_list(agents), to_list(products)


def _analyze_smiles(smiles: str) -> Dict:
    _require_rdkit()
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem import rdMolDescriptors

    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return {"input_smiles": smiles, "ok": False, "error": "RDKit failed to parse"}

    formula = rdMolDescriptors.CalcMolFormula(mol) or ""
    mw = Descriptors.MolWt(mol)
    heavy = mol.GetNumHeavyAtoms()
    rings = rdMolDescriptors.CalcNumRings(mol)
    arom_rings = rdMolDescriptors.CalcNumAromaticRings(mol)

    elem_counts: Dict[str, int] = {}
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        elem_counts[sym] = elem_counts.get(sym, 0) + 1

    has_map = any(a.GetAtomMapNum() > 0 for a in mol.GetAtoms())

    return {
        "input_smiles": smiles,
        "ok": True,
        "formula": formula,
        "mw": float(mw),
        "heavy_atoms": int(heavy),
        "rings": int(rings),
        "aromatic_rings": int(arom_rings),
        "element_counts": elem_counts,
        "has_atom_map": bool(has_map),
    }


def _merge_counts(items: List[Dict]) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for it in items:
        if it.get("ok"):
            for k, v in it.get("element_counts", {}).items():
                counts[k] = counts.get(k, 0) + int(v)
    return counts


def analyze_reaction_smiles(rxn_smiles: str) -> Dict:
    """
    Basic analysis of a reaction SMILES string.
    Input format: reactants>agents>products (agents optional)
    Returns a dict summarizing per-molecule info and basic element counts.
    """
    reactants, agents, products = _split_reaction_smiles(rxn_smiles)

    r_infos = [_analyze_smiles(s) for s in reactants]
    a_infos = [_analyze_smiles(s) for s in agents]
    p_infos = [_analyze_smiles(s) for s in products]

    atom_map_present = any(
        info.get("has_atom_map") for info in (r_infos + a_infos + p_infos) if info.get("ok")
    )

    result = {
        "reactants": r_infos,
        "agents": a_infos,
        "products": p_infos,
        "atom_map_present": bool(atom_map_present),
        "element_counts": {
            "reactants": _merge_counts(r_infos),
            "products": _merge_counts(p_infos),
        },
    }
    return result

