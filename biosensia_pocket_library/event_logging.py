"""Structured, scrubbed dual-sink run logging."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .scrub import scrub


class EventLogger:
    def __init__(self, run_dir: Path):
        log_dir = run_dir / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        self.text_path = log_dir / "pipeline.log"
        self.json_path = log_dir / "events.jsonl"

    def emit(self, level: str, stage: str, code: str, message: str, *, complex_id: str | None = None,
             pocket_instance_id: str | None = None, details: dict[str, Any] | None = None) -> None:
        event = scrub({"timestamp_utc": datetime.now(timezone.utc).isoformat(), "level": level,
                       "stage": stage, "complex_id": complex_id, "pocket_instance_id": pocket_instance_id,
                       "event_code": code, "message": message, "worker_id": 0, "details": details or {}})
        with self.json_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(event, sort_keys=True, ensure_ascii=False) + "\n")
        with self.text_path.open("a", encoding="utf-8") as handle:
            handle.write(f"{event['timestamp_utc']} {level.upper()} {stage} {code} {event['message']}\n")
