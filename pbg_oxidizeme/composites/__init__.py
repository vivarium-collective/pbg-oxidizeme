"""Composite generators for pbg-oxidizeme."""

from . import steady_state  # noqa: F401 — register decorator side effect

from .steady_state import oxidizeme_steady_state

__all__ = ["oxidizeme_steady_state"]
