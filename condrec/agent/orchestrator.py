from __future__ import annotations

import json
import time
from typing import Dict, Generator, List, Optional

from ..analysis import analyze_reaction_smiles
from ..config import get_settings
from ..llm import LLMClient, build_initial_messages
from ..search import build_index, search_similar
from ..schemas import (
    LLMRecommendation,
    OrchestratorEvent,
    OrchestratorResult,
    ReactionAnalysis,
    ReactionInput,
)
from ..validation.checks import validate_recommendations


class Orchestrator:
    """Coordinates local tools with the LLM in a loop.

    The primary APIs are:
      - iter_steps(): yields events for streaming to a UI
      - run_once(): performs a single recommendation pass and returns a result
    """

    def __init__(self):
        self.settings = get_settings()
        # Build similarity index lazily at first use
        self._index_built = False
        self._index = None
        self._llm = LLMClient()

    def _ensure_index(self):
        if not self._index_built:
            self._index = build_index(self.settings.dataset_path)
            self._index_built = True

    def iter_steps(
        self, rxn: ReactionInput, *, top_k: int = 5
    ) -> Generator[OrchestratorEvent, None, OrchestratorResult]:
        t0 = time.time()
        yield {"type": "status", "msg": "Analyzing reaction"}
        analysis: ReactionAnalysis = analyze_reaction_smiles(rxn.raw_smiles)
        yield {"type": "analysis", "data": analysis}

        yield {"type": "status", "msg": "Searching local precedents"}
        self._ensure_index()
        similar = search_similar(self._index, rxn.raw_smiles, k=top_k) if self._index else []
        yield {"type": "similar", "data": similar}

        summary = {
            "analysis": analysis,
            "input": {"smiles": rxn.raw_smiles, **({} if not rxn.meta else {"meta": rxn.meta})},
        }
        messages = build_initial_messages(summary, similar)
        yield {"type": "llm_request", "data": messages}

        resp = self._llm.generate(messages)
        yield {"type": "llm_response", "data": resp}

        # Try to parse JSON out of the model's response
        text = resp.get("text") or ""
        try:
            obj = json.loads(text)
        except Exception:
            obj = {"conditions": []}

        recs: List[LLMRecommendation] = list(obj.get("conditions") or [])
        issues = validate_recommendations(recs, analysis)
        yield {"type": "validation", "data": issues}

        elapsed_ms = int((time.time() - t0) * 1000)
        yield {
            "type": "recommendations",
            "data": {
                "analysis": analysis,
                "similar": similar,
                "recommendations": recs,
                "issues": issues,
                "elapsed_ms": elapsed_ms,
            },
        }

        return {
            "analysis": analysis,
            "similar_hits": similar,
            "llm_payload": {"reaction_summary": summary, "similar_hits": similar},
            "raw_llm_output": resp,
            "recommendations": recs,
            "validation_issues": issues,
            "elapsed_ms": elapsed_ms,
        }

    def run_once(self, rxn: ReactionInput, *, top_k: int = 5) -> OrchestratorResult:
        gen = self.iter_steps(rxn, top_k=top_k)
        last: Optional[OrchestratorResult] = None
        for _ in gen:
            # drain events to completion
            pass
        try:
            last = gen.send(None)  # get the return value from generator
        except StopIteration as e:  # type: ignore
            last = e.value  # type: ignore
        return last  # type: ignore
