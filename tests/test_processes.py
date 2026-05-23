"""Shape tests for OxidizeMEStep.

The real-bridge solve is gated behind upstream packages (cobrame, ecolime,
qminospy, oxidizeme) — see the README. These tests check the wrapper's
interface contract: ports, config, error semantics. The end-to-end solve
test ``test_real_bridge_solve`` uses ``pytest.importorskip`` so CI without
the upstream stack stays green; on a dev machine with the stack installed,
it actually exercises the real bridge.
"""

from __future__ import annotations

import pytest


def test_step_class_is_importable():
    from pbg_oxidizeme import OxidizeMEStep, OxidizeMEMissingError

    assert OxidizeMEStep is not None
    assert issubclass(OxidizeMEMissingError, RuntimeError)


def test_step_ports_shape():
    from pbg_oxidizeme import OxidizeMEStep
    from process_bigraph import allocate_core

    core = allocate_core()
    step = OxidizeMEStep(config={}, core=core)
    inputs = step.inputs()
    outputs = step.outputs()

    # Cell-side interface contract surface
    assert set(inputs) == {
        "external_concentrations",
        "ros_concentrations",
        "metal_concentrations",
        "mu_fixed",
        "atpm_mmol_per_gDW_per_h",
    }
    assert set(outputs) == {
        "mu",
        "external_exchange_fluxes",
        "proteome_allocations",
        "damaged_proteome_fraction",
        "solver_status",
        "wall_time_s",
    }


def test_initial_state_carries_defaults():
    from pbg_oxidizeme import OxidizeMEStep
    from process_bigraph import allocate_core

    core = allocate_core()
    step = OxidizeMEStep(config={}, core=core)
    init = step.initial_state()

    assert "ros_concentrations" in init
    assert "h2o2_c" in init["ros_concentrations"]
    assert init["ros_concentrations"]["h2o2_c"] > 0


def test_strict_raises_when_upstream_missing():
    """When the OxidizeME stack isn't installed, strict mode (default) must
    raise OxidizeMEMissingError with the install pointer in the message.
    """
    from pbg_oxidizeme import OxidizeMEStep, OxidizeMEMissingError
    from process_bigraph import allocate_core

    core = allocate_core()
    step = OxidizeMEStep(config={"strict": True}, core=core)
    state = step.initial_state()
    try:
        import oxidizeme  # noqa: F401
    except ImportError:
        with pytest.raises(OxidizeMEMissingError) as ei:
            step.update(state)
        msg = str(ei.value)
        assert "OxidizeME" in msg or "oxidizeme" in msg
        assert "install" in msg.lower() or "Docker" in msg or "qMINOS" in msg
    else:
        pytest.skip("OxidizeME is installed locally — strict-error path not exercised.")


def test_non_strict_returns_labeled_empty_when_upstream_missing():
    """Non-strict mode (dashboard-scaffolding path) must NOT raise; it
    returns a labelled empty solution with ``solver_status='upstream_missing'``.
    """
    from pbg_oxidizeme import OxidizeMEStep
    from process_bigraph import allocate_core

    core = allocate_core()
    step = OxidizeMEStep(config={"strict": False}, core=core)
    state = step.initial_state()
    try:
        import oxidizeme  # noqa: F401
        pytest.skip("OxidizeME is installed locally — upstream-missing path not exercised.")
    except ImportError:
        out = step.update(state)
        assert out["solver_status"] == "upstream_missing"
        assert out["external_exchange_fluxes"] == {}
        # mu / damaged_proteome are NaN (we don't fake science)
        import math
        assert math.isnan(out["mu"])
        assert math.isnan(out["damaged_proteome_fraction"])


def test_real_bridge_solve():
    """End-to-end real-bridge test. Skips when the upstream stack is absent."""
    oxidizeme = pytest.importorskip("oxidizeme")
    qminospy  = pytest.importorskip("qminospy")
    # If we get here, the upstream stack is installed. We still need a
    # pickled model — without one the test is moot.
    import os
    pkl = os.environ.get("OXIDIZEME_PICKLE")
    if not pkl or not os.path.exists(pkl):
        pytest.skip("Set OXIDIZEME_PICKLE to a pickled StressME model to exercise the real bridge.")

    from pbg_oxidizeme import OxidizeMEStep
    from process_bigraph import allocate_core

    core = allocate_core()
    step = OxidizeMEStep(
        config={
            "me_model_pickle_path": pkl,
            "strict": True,
        },
        core=core,
    )
    state = step.initial_state()
    state["external_concentrations"] = {
        "GLC[p]":             22.2,
        "OXYGEN-MOLECULE[p]": 0.21,
        "AMMONIUM[p]":        19.0,
    }
    out = step.update(state)
    # Real solve should produce a finite μ near literature 0.6–0.7 1/h
    assert out["solver_status"] in {"optimal", "infeasible"}
    if out["solver_status"] == "optimal":
        assert 0.05 < out["mu"] < 2.0
        assert isinstance(out["external_exchange_fluxes"], dict)
