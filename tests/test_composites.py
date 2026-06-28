"""Composite-generator tests for pbg-oxidizeme.

Asserts that the steady-state generator registers via the
@composite_generator decorator and that the generated document has the
expected store layout — without requiring the upstream OxidizeME stack
to run (the generator is just dict-building).
"""

from __future__ import annotations

import pytest


def test_generator_is_registered():
    # Importing the package's composites subpackage runs the @composite_generator
    # decorator side effect — this is the documented "import to register" pattern.
    import pbg_oxidizeme.composites  # noqa: F401
    from pbg_superpowers.composite_generator import _REGISTRY

    matches = [eid for eid in _REGISTRY if eid.endswith(".oxidizeme_steady_state")]
    assert matches, f"oxidizeme_steady_state missing; have {sorted(_REGISTRY)[:8]}"


def test_generated_document_has_expected_stores():
    from pbg_oxidizeme.composites import oxidizeme_steady_state

    doc = oxidizeme_steady_state(
        core=None,
        me_model_pickle_path="",
        glucose_mM=22.2,
        oxygen_mM=0.21,
        ammonium_mM=19.0,
    )
    assert "oxidizeme" in doc
    assert "stores" in doc
    assert "emitter" in doc

    stores = doc["stores"]
    assert "external_concentrations" in stores
    assert "ros_concentrations" in stores
    assert stores["external_concentrations"]["GLC[p]"] == 22.2
    assert stores["external_concentrations"]["OXYGEN-MOLECULE[p]"] == 0.21


def test_composite_assembles_in_non_strict_mode():
    """Build the Composite without running it. Exercises schema reconciliation
    end-to-end on the generated document. Non-strict so the wrapper doesn't
    need OxidizeME installed.
    """
    from pbg_oxidizeme import OxidizeMEStep
    from pbg_oxidizeme.composites import oxidizeme_steady_state
    from process_bigraph import Composite, allocate_core

    core = allocate_core()
    # The generated document wires the process by its `local:OxidizeMEStep`
    # address, resolved from the core's link registry. In a dashboard the
    # workspace's processes are discovered/registered for us; in an isolated
    # test venv they are not, so register defensively (idempotent — a dict
    # assignment) before assembling the Composite.
    core.register_link("OxidizeMEStep", OxidizeMEStep)
    doc = oxidizeme_steady_state(core=core, strict=False)
    sim = Composite({"state": doc}, core=core)
    assert sim is not None
