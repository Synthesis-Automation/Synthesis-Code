from __future__ import annotations

import json
from typing import Any, Dict, List, Optional

from ..config import get_settings


class LLMClient:
    """Pluggable LLM client with a development-friendly dry-run mode.

    - When CONDREC_LLM_DRYRUN is enabled (default), returns a canned response.
    - Otherwise, attempts to call the configured provider.
    """

    def __init__(self, provider: Optional[str] = None, model: Optional[str] = None):
        s = get_settings()
        self.provider = provider or s.provider
        self.model = model or s.model
        self.temperature = s.temperature
        self.top_p = s.top_p
        self.seed = s.seed
        self._dry_run = s.dry_run

    def generate(self, messages: List[Dict[str, str]], **kwargs) -> Dict[str, Any]:
        if self._dry_run:
            return self._fake_response()
        if self.provider == "openai":
            return self._call_openai(messages, **kwargs)
        elif self.provider == "anthropic":
            return self._call_anthropic(messages, **kwargs)
        else:
            raise RuntimeError(f"Unsupported provider: {self.provider}")

    # --- Providers --------------------------------------------------------

    def _call_openai(self, messages: List[Dict[str, str]], **_: Any) -> Dict[str, Any]:  # pragma: no cover
        try:
            import openai  # type: ignore
        except Exception as e:
            raise RuntimeError(
                "OpenAI client not installed. 'pip install openai' or set CONDREC_LLM_DRYRUN=1"
            ) from e

        s = get_settings()
        if s.openai_base_url:
            openai.base_url = s.openai_base_url
        if s.openai_api_key:
            openai.api_key = s.openai_api_key
        else:
            raise RuntimeError("OPENAI_API_KEY is required or enable CONDREC_LLM_DRYRUN=1")

        resp = openai.chat.completions.create(
            model=self.model,
            messages=messages,
            temperature=self.temperature,
            top_p=self.top_p,
        )
        # Normalize
        content = resp.choices[0].message.content if resp.choices else ""
        return {"text": content, "raw": resp.model_dump() if hasattr(resp, "model_dump") else resp}

    def _call_anthropic(self, messages: List[Dict[str, str]], **_: Any) -> Dict[str, Any]:  # pragma: no cover
        try:
            import anthropic  # type: ignore
        except Exception as e:
            raise RuntimeError(
                "Anthropic client not installed. 'pip install anthropic' or set CONDREC_LLM_DRYRUN=1"
            ) from e

        s = get_settings()
        if not s.anthropic_api_key:
            raise RuntimeError("ANTHROPIC_API_KEY is required or enable CONDREC_LLM_DRYRUN=1")

        client = anthropic.Anthropic(api_key=s.anthropic_api_key)
        # Convert OpenAI-style roles to Anthropic messages
        sys_prompt = ""
        user_content = []
        for m in messages:
            if m.get("role") == "system":
                sys_prompt = m.get("content", "")
            elif m.get("role") == "user":
                user_content.append({"type": "text", "text": m.get("content", "")})

        resp = client.messages.create(
            model=self.model,
            max_tokens=1000,
            temperature=self.temperature,
            system=sys_prompt or None,
            messages=[{"role": "user", "content": user_content}],
        )
        text = "".join(part.text for part in resp.content if getattr(part, "type", "text") == "text")
        return {"text": text, "raw": getattr(resp, "model_dump", lambda: resp)()}

    # --- Helpers ---------------------------------------------------------

    def _fake_response(self) -> Dict[str, Any]:
        """Return deterministic, plausible JSON for development."""
        fake = {
            "conditions": [
                {
                    "solvent": "THF",
                    "reagent": "NaH",
                    "base": None,
                    "catalyst": None,
                    "temperature_c": 25.0,
                    "time_h": 2.0,
                    "atmosphere": "N2",
                    "concentration_m": 0.2,
                    "stoichiometry": {"reagent": 1.2},
                    "steps": ["Charge reagents", "Stir at RT for 2 h"],
                    "rationale": "Baseline conditions; mild base in aprotic solvent.",
                    "references": [],
                },
                {
                    "solvent": "DMF",
                    "reagent": "K2CO3",
                    "temperature_c": 80.0,
                    "time_h": 4.0,
                    "atmosphere": "Ar",
                    "stoichiometry": {"reagent": 2.0},
                    "steps": ["Heat to 80C"],
                    "rationale": "Hot DMF variant to drive conversion.",
                    "references": [],
                },
            ]
        }
        return {"text": json.dumps(fake), "raw": fake}

