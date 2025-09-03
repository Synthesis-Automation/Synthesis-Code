from __future__ import annotations

from typing import Any, Dict, List


def build_initial_messages(summary: Dict[str, Any], similar_hits: List[Dict[str, Any]]) -> List[Dict[str, str]]:
    """Construct a chat-style message list for a general-purpose LLM.

    The content is provider-agnostic (role/content pairs).
    """
    sys_msg = (
        "You are an expert synthetic chemist assisting with reaction condition "
        "recommendations. Propose concrete, plausible conditions with brief rationale."
    )

    user_msg = (
        "Analyze the following reaction and propose 3 candidate condition sets.\n\n"
        "Reaction summary (JSON):\n"
        f"{summary}\n\n"
        "Top similar precedents (JSON):\n"
        f"{similar_hits}\n\n"
        "Output a concise JSON with fields: conditions[], each with solvent, reagent/catalyst/base, "
        "temperature_c, time_h, atmosphere (if relevant), and a short rationale and references."
    )

    return [
        {"role": "system", "content": sys_msg},
        {"role": "user", "content": user_msg},
    ]

