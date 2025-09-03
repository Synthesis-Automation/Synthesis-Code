from __future__ import annotations

from typing import List

from ..schemas import LLMRecommendation, ReactionAnalysis


def _temperature_ok(t: float) -> bool:
    return -100.0 <= t <= 250.0


def _time_ok(h: float) -> bool:
    return 0.0 < h <= 168.0


def validate_recommendations(recs: List[LLMRecommendation], analysis: ReactionAnalysis) -> List[str]:
    """Shallow validation checks for obviously invalid suggestions.

    Returns a list of human-readable issues (empty if none).
    """
    issues: List[str] = []
    for i, r in enumerate(recs, start=1):
        t = r.get("temperature_c")
        if t is not None and not _temperature_ok(float(t)):
            issues.append(f"#{i}: temperature_c {t}C out of typical bounds")
        h = r.get("time_h")
        if h is not None and not _time_ok(float(h)):
            issues.append(f"#{i}: time_h {h}h out of typical bounds")
        conc = r.get("concentration_m")
        if conc is not None and (float(conc) <= 0 or float(conc) > 5):
            issues.append(f"#{i}: concentration {conc}M looks unreasonable")
    return issues
