"""Central redaction helpers used before persistence or logging."""

from __future__ import annotations

import re
from typing import Any

PDF_TOKEN = re.compile(r"(?i)(?<!\S)[^\s]*?\.pdf(?:[\]\[(){}.,;:]*)?(?!\S)")


def remove_pdf_tokens(text: str) -> str:
    """Remove whitespace-delimited tokens whose path basename ends in .pdf."""
    kept: list[str] = []
    remove_next_separator = False
    for token in text.split():
        if remove_next_separator and re.fullmatch(r"[,;:|/\\-]+", token):
            remove_next_separator = False
            continue
        remove_next_separator = False
        core = token.strip("()[]{}<>,;:'\"")
        basename = re.split(r"[/\\]", core)[-1]
        if basename.lower().endswith(".pdf"):
            if kept and re.fullmatch(r"[,;:|/\\-]+", kept[-1]):
                kept.pop()
            remove_next_separator = True
            continue
        kept.append(token)
    result = " ".join(kept)
    result = re.sub(r"\s+([,;:])", r"\1", result)
    return result.strip(" \t,;:")


def scrub(value: Any) -> Any:
    if isinstance(value, str):
        return remove_pdf_tokens(value)
    if isinstance(value, dict):
        return {str(key): scrub(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [scrub(item) for item in value]
    return value
