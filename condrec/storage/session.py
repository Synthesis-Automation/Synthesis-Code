from __future__ import annotations

import json
import os
import time
import uuid
from typing import Any, Dict

from ..config import get_settings


class SessionLogger:
    """Append-only JSONL session logger for debugging/evaluation."""

    def __init__(self, session_id: str | None = None):
        self.settings = get_settings()
        self.session_id = session_id or time.strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:8]
        self.path = os.path.join(self.settings.session_dir, f"{self.session_id}.jsonl")

    def log(self, event: Dict[str, Any]) -> None:
        try:
            with open(self.path, "a", encoding="utf-8") as f:
                f.write(json.dumps(event, ensure_ascii=False) + "\n")
        except Exception:
            pass

