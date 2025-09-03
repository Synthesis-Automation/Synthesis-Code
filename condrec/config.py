"""Configuration helpers for CondRec.

Environment-first with light defaults. Avoids heavy config systems.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional


def _getenv_float(name: str, default: float) -> float:
    try:
        return float(os.getenv(name, str(default)))
    except Exception:
        return default


def _getenv_int(name: str, default: int) -> int:
    try:
        return int(os.getenv(name, str(default)))
    except Exception:
        return default


@dataclass
class Settings:
    # Data / cache
    dataset_path: Optional[str] = os.getenv("CONDREC_DATASET_PATH")
    cache_dir: str = os.getenv("CONDREC_CACHE_DIR", ".condrec_cache")
    session_dir: str = os.getenv("CONDREC_SESSION_DIR", ".condrec_sessions")

    # LLM provider & model
    provider: str = os.getenv("CONDREC_LLM_PROVIDER", "openai")  # openai|anthropic|local
    model: str = os.getenv("CONDREC_LLM_MODEL", "gpt-4o-mini")
    temperature: float = _getenv_float("CONDREC_TEMPERATURE", 0.2)
    top_p: float = _getenv_float("CONDREC_TOP_P", 1.0)
    seed: Optional[int] = None if os.getenv("CONDREC_SEED") is None else _getenv_int("CONDREC_SEED", 0)

    # Provider creds / endpoints
    openai_api_key: Optional[str] = os.getenv("OPENAI_API_KEY")
    openai_base_url: Optional[str] = os.getenv("OPENAI_BASE_URL")
    anthropic_api_key: Optional[str] = os.getenv("ANTHROPIC_API_KEY")

    # Dry-run mode to avoid network calls during development
    dry_run: bool = os.getenv("CONDREC_LLM_DRYRUN", "1") not in ("0", "false", "False")


_SETTINGS: Optional[Settings] = None


def get_settings() -> Settings:
    global _SETTINGS
    if _SETTINGS is None:
        _SETTINGS = Settings()
        # Ensure directories exist
        for path in (_SETTINGS.cache_dir, _SETTINGS.session_dir):
            try:
                os.makedirs(path, exist_ok=True)
            except Exception:
                pass
    return _SETTINGS

