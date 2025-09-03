"""Core domain types for the conditions recommendation workflow.

These types are intentionally light-weight (dataclasses and TypedDicts)
so the rest of the app can import them without heavy dependencies.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Literal, Optional, TypedDict


class MoleculeInfo(TypedDict, total=False):
    input_smiles: str
    ok: bool
    formula: str
    mw: float
    heavy_atoms: int
    rings: int
    aromatic_rings: int
    element_counts: Dict[str, int]
    has_atom_map: bool
    error: str


class ReactionAnalysis(TypedDict):
    reactants: List[MoleculeInfo]
    agents: List[MoleculeInfo]
    products: List[MoleculeInfo]
    atom_map_present: bool
    # {"reactants": {"C": 10, ...}, "products": {...}}
    element_counts: Dict[str, Dict[str, int]]


@dataclass(frozen=True)
class ReactionInput:
    """Minimal representation of a reaction request from the user/UI."""

    raw_smiles: str
    # Optional user-specified metadata (e.g., name, project, notes)
    meta: Dict[str, Any] = field(default_factory=dict)


class SimilarityHit(TypedDict, total=False):
    source_id: str
    reaction_smiles: str
    score: float
    metadata: Dict[str, Any]


class LLMQueryPayload(TypedDict, total=False):
    reaction_summary: Dict[str, Any]
    similar_hits: List[SimilarityHit]
    constraints: Dict[str, Any]
    user_feedback: Optional[str]
    output_schema: Dict[str, Any]


class LLMRecommendation(TypedDict, total=False):
    # Generic condition slots – extend per domain as needed
    solvent: Optional[str]
    catalyst: Optional[str]
    base: Optional[str]
    reagent: Optional[str]
    temperature_c: Optional[float]
    time_h: Optional[float]
    atmosphere: Optional[str]
    concentration_m: Optional[float]
    stoichiometry: Dict[str, float]
    steps: List[str]
    rationale: str
    references: List[str]


class OrchestratorEvent(TypedDict, total=False):
    type: Literal[
        "status",
        "analysis",
        "similar",
        "llm_request",
        "llm_response",
        "recommendations",
        "validation",
        "error",
    ]
    data: Any
    msg: Optional[str]


class OrchestratorResult(TypedDict, total=False):
    analysis: ReactionAnalysis
    similar_hits: List[SimilarityHit]
    llm_payload: LLMQueryPayload
    raw_llm_output: Dict[str, Any]
    recommendations: List[LLMRecommendation]
    validation_issues: List[str]
    elapsed_ms: int

