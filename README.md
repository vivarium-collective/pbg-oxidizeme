# pbg-oxidizeme

A process-bigraph wrapper for **OxidizeME** — a genome-scale model of
*Escherichia coli* metabolism + macromolecular expression (ME-model) with
reactive-oxygen-species damage and repair reactions, published by Yang et al.
([2019 PNAS](https://doi.org/10.1073/pnas.1905039116), 2018 bioRxiv preprint)
from the [Palsson lab](https://systemsbiology.ucsd.edu) at UCSD's Systems
Biology Research Group.

The upstream OxidizeME model is the engine used in Beulig et al. (2025,
mSystems, [doi:10.1128/msystems.00323-25](https://doi.org/10.1128/msystems.00323-25))
to interpret high-cell-density fed-batch *E. coli* physiology. This wrapper
exposes OxidizeME as a process-bigraph `Step`, conforming to the
**cell-side interface contract** of the [v2ecoli
multiscale-bioprocess](https://github.com/vivarium-collective/v2ecoli/tree/multiscale-bioprocess)
investigation, so it can serve as a candidate-future engine alongside
v2ecoli baseline.

## What it does

`OxidizeMEStep` takes the cell-side-interface inputs:

- `external_concentrations` (`map[float]`, mM) — environment (glucose, O2, …)
- `ros_concentrations` (`map[float]`, nmol/L) — intracellular H₂O₂, O₂⁻⁻
- `metal_concentrations` (`map[float]`, nmol/L) — Fe²⁺, Mn²⁺, Zn²⁺
- `mu_fixed` (`maybe[float]`, 1/h) — if set, fixed-μ LP; else bisection for μ_max
- `atpm_mmol_per_gDW_per_h` (`maybe[float]`) — ATP maintenance bound override

…and emits the optimization outputs:

- `mu` (`overwrite[float]`, 1/h)
- `external_exchange_fluxes` (`map[float]`, mmol/(gDW·h))
- `proteome_allocations` (`map[float]`, fraction per pathway)
- `damaged_proteome_fraction` (`overwrite[float]`, 0..1)
- `solver_status` (`overwrite[string]`)
- `wall_time_s` (`overwrite[float]`)

Each call solves the ME-LP once — bisection for μ_max by default, fixed-μ LP
when `mu_fixed` is provided. The Step is the natural shape because ME-models
are constraint-based optimizations, not time-stepped simulations. A surrounding
bigraph that wants time dynamics wires this Step downstream of a Process
that holds the time integration (the dFBA pattern; see the v2ecoli mbp-03
study for the worked example).

## Real bridge — install obstacle

This package is a **real bridge**, not a reproduction or mock. It drives the
genuine OxidizeME / cobrame / ecolime / solvemepy / qminospy stack. That stack
is **environmentally constrained**:

1. **Python 2.7 / cobrapy 0.5.11 era.** Upstream cobrame, ecolime, and
   oxidizeme were last touched in 2018–2019. A Python-3 port exists on the
   `devel` branches but is not on PyPI.
2. **qMINOS — proprietary 128-bit-precision solver.** Obtained by contacting
   **Prof. M. A. Saunders, Stanford University**. ME-models are inherently
   ill-scaled; only extended-precision solvers (qMINOS, or the 80-bit
   academic-free SoPlex) converge reliably.
3. **The ME-model itself is heavy.** Build-from-scratch via ecolime takes
   minutes; the recommended path is to load a pickled `StressME` model.

### Recommended install path

The SBRG team publishes Docker images at
[`hub.docker.com/r/sbrg/cobrame`](https://hub.docker.com/r/sbrg/cobrame) with
the full stack pre-installed. Build a wrapper environment around the Docker
runtime (e.g. via `docker exec`) and supply a path to a pickled `StressME`
model:

```python
from pbg_oxidizeme import OxidizeMEStep
from process_bigraph import allocate_core

core = allocate_core()
step = OxidizeMEStep(
    config={
        "me_model_pickle_path": "/path/to/stress_me.pickle",
        "strict": True,
    },
    core=core,
)
state = step.initial_state()
state["external_concentrations"] = {
    "GLC[p]":             22.2,    # mM, M9 ~4 g/L
    "OXYGEN-MOLECULE[p]": 0.21,    # mM, near-saturation aerobic
    "AMMONIUM[p]":        19.0,    # mM, M9 ~1 g/L NH4Cl
}
solution = step.update(state)
print(solution["mu"], "1/h")
```

### Dashboard-scaffolding path (no solver)

For wiring composites in the dashboard before the solver is set up, instantiate
in non-strict mode:

```python
step = OxidizeMEStep(config={"strict": False}, core=core)
out = step.update(state)
# out["solver_status"] == "upstream_missing"
# out["mu"] is NaN — we don't fake science
```

This returns a labelled empty solution (`solver_status: "upstream_missing"`,
NaN scalars, empty dicts) so the dashboard's Composites tab can browse and
render the Step's port surface without the solver. **Never assert numerical
fidelity on those outputs** — they are deliberately inert.

## Installation

```bash
# From source (development):
git clone https://github.com/vivarium-collective/pbg-oxidizeme
cd pbg-oxidizeme
uv venv .venv && source .venv/bin/activate
uv pip install -e ".[dev]"
pytest
```

The `bigraph-schema` + `process-bigraph` deps install at this step. The
OxidizeME stack (cobrame, ecolime, qminospy, oxidizeme) is **NOT** a pip
dep — see "Real bridge — install obstacle" above for that.

> Once installed, processes register automatically via
> `bigraph_schema.package.discover` — no manual `register_link()` calls
> are needed.

## Demo

```bash
python demo/demo_report.py
```

The demo runs the composite generator in non-strict mode (so it works
without the real solver) and produces `demo/report.html` with the Step's
port surface, an architecture diagram, and the labelled empty solution.
When OxidizeME is installed, the same demo exercises the real bridge.

## Cell-side interface contract

`OxidizeMEStep` is designed to satisfy the substitutability contract documented
at
[v2ecoli/.pbg/worktrees/multiscale-bioprocess/references/expert/cell_side_interface_contract.md](https://github.com/vivarium-collective/v2ecoli/blob/multiscale-bioprocess/references/expert/cell_side_interface_contract.md).
In particular:

- Input port names + types match the contract surface
- External-molecule keys use EcoCyc-style `<MOL>[p]` IDs (mapped to the
  ME-model's BiGG IDs via `config.exchange_id_map`)
- The Step is drop-in interchangeable with `v2ecoli.composites.baseline`
  under the same reactor composite — wired into the same `environment.*`
  + `agents.*` store paths

When mbp-06's gap-analysis recommends a v2ecoli-vs-OxidizeME comparator
study, this wrapper is the engine half of that comparison.

## Limitations

- **Bound-mapping calibration.** The exchange-bound mapping in
  `_apply_exchange_bounds` is a placeholder Monod-style soft cap with
  literature-anchored v_max defaults. The v2ecoli mbp investigation's
  spec PR replaces it with a per-substrate calibrated form.
- **Proteome-allocation classifier.** `_compute_proteome_allocations`
  groups translation fluxes into glycolysis / PPP / TCA / other by a
  keyword match on linked reaction IDs. This is approximate; the Beulig
  2025 analysis uses curated pathway membership and iModulon assignments.
  Override via a future config field.
- **No fed-batch operations.** OxidizeME is a steady-state ME-model; the
  fed-batch surrounding logic (feed schedules, volume dynamics) lives in
  the encapsulating composite, not this Step.
- **Strain scope.** Defaults target WT *E. coli* K-12 MG1655 (the iJL1678b-ME
  base). SGKO / TRP / MEL strain comparisons (from Beulig 2025) require
  separate pickled models with the appropriate gene knockouts /
  heterologous-pathway expansions wired in.

## Attribution

OxidizeME (the upstream model) is © Laurence Yang and the Systems Biology
Research Group, UCSD. This wrapper does not redistribute the OxidizeME
source — it provides a process-bigraph adapter that drives a locally-
installed copy. Please cite Yang et al. when publishing results derived
from this wrapper:

> Yang L, Mih N, Anand A, Park JH, Tan J, Yurkovich JT, Monk JM, Lloyd CJ,
> Sandberg TE, Seo SW, Kim D, Sastry AV, Phaneuf P, Gao Y, Broddrick JT,
> Chen K, Heckmann D, Szubin R, Hefner Y, Feist AM, Palsson BO. 2019.
> Cellular responses to reactive oxygen species are predicted from
> molecular mechanisms. *Proc Natl Acad Sci USA* **116**:14368–14373.

## License

MIT. See `LICENSE`.
