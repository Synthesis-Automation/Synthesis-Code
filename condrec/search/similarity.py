from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

from ..schemas import SimilarityHit


def _require_rdkit():
    try:
        from rdkit import Chem  # noqa: F401
        from rdkit.Chem import rdMolDescriptors, DataStructs  # noqa: F401
    except Exception as e:  # pragma: no cover - runtime dependency
        raise ImportError(
            "RDKit is required for similarity search.\n"
            "Install via: conda install -c conda-forge rdkit"
        ) from e


def _split_reaction_smiles(rxn: str) -> Tuple[List[str], List[str], List[str]]:
    parts = rxn.strip().split(">")
    if len(parts) == 2:
        parts = [parts[0], "", parts[1]]
    if len(parts) != 3:
        raise ValueError("Reaction SMILES must have 2 or 3 parts separated by '>'")
    to_list = lambda s: [p for p in s.split(".") if p] if s else []
    return to_list(parts[0]), to_list(parts[1]), to_list(parts[2])


def _morgan_fp(smiles: str, radius: int = 2, n_bits: int = 2048):
    _require_rdkit()
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return None
    return rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)


def _rxn_fp(rxn_smiles: str, radius: int = 2, n_bits: int = 2048):
    """Simple reaction-level fingerprint: XOR of reactant/product Morgan FPs.

    This is a pragmatic placeholder until a dedicated reaction FP is selected.
    """
    _require_rdkit()
    from rdkit.Chem import DataStructs

    r, a, p = _split_reaction_smiles(rxn_smiles)
    fps = []
    for s in (r + p):  # ignore agents for similarity baseline
        fp = _morgan_fp(s, radius=radius, n_bits=n_bits)
        if fp is not None:
            fps.append(fp)
    if not fps:
        return None
    acc = fps[0]
    for fp in fps[1:]:
        acc ^= fp  # bitwise XOR to capture change
    return acc


@dataclass
class Record:
    source_id: str
    reaction_smiles: str
    metadata: Dict


class SimilarityIndex:
    def __init__(self, radius: int = 2, n_bits: int = 2048):
        self.radius = radius
        self.n_bits = n_bits
        self._records: List[Record] = []
        self._fps = []  # RDKit ExplicitBitVect list

    def add(self, rec: Record):
        fp = _rxn_fp(rec.reaction_smiles, self.radius, self.n_bits)
        if fp is None:
            return
        self._records.append(rec)
        self._fps.append(fp)

    def add_many(self, records: Iterable[Record]):
        for r in records:
            self.add(r)

    def search(self, rxn_smiles: str, k: int = 10) -> List[SimilarityHit]:
        _require_rdkit()
        from rdkit.Chem import DataStructs

        qfp = _rxn_fp(rxn_smiles, self.radius, self.n_bits)
        if qfp is None or not self._fps:
            return []
        scores = [DataStructs.TanimotoSimilarity(qfp, fp) for fp in self._fps]
        idxs = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)[:k]
        hits: List[SimilarityHit] = []
        for i in idxs:
            rec = self._records[i]
            hits.append(
                {
                    "source_id": rec.source_id,
                    "reaction_smiles": rec.reaction_smiles,
                    "score": float(scores[i]),
                    "metadata": rec.metadata,
                }
            )
        return hits


def build_index(dataset_path: Optional[str]) -> SimilarityIndex:
    """Build an index from a simple JSONL dataset if available.

    Expected JSONL schema per line (flexible):
      {"id": "...", "reaction_smiles": "...", "metadata": {...}}
    If path is missing or file unreadable, returns an empty index.
    """
    import json
    import os

    index = SimilarityIndex()
    if not dataset_path:
        return index
    if not os.path.exists(dataset_path):
        return index

    try:
        with open(dataset_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                obj = json.loads(line)
                rid = str(obj.get("id") or obj.get("uid") or obj.get("_id") or len(index._records))
                rxn = obj.get("reaction_smiles") or obj.get("rxn") or obj.get("reaction")
                if not rxn:
                    continue
                meta = obj.get("metadata") or {k: v for k, v in obj.items() if k not in ("id", "reaction_smiles")}
                index.add(Record(rid, rxn, meta))
    except Exception:
        # Be tolerant; an empty index is acceptable in early development.
        return index

    return index


def search_similar(index: SimilarityIndex, rxn_smiles: str, k: int = 10) -> List[SimilarityHit]:
    return index.search(rxn_smiles, k=k)
