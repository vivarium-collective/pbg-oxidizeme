"""OxidizeMEStep — process-bigraph Step that drives the real OxidizeME ME-model.

The wrapped tool is StressME (oxidizeme.model.StressME), which takes an
ME-model NLP wrapper (qminospy.me1.ME_NLP1 around a cobrame.MEModel built
via ecolime.build_ME_model) and adds reactive-oxygen-species (ROS) damage
and repair reactions on top. A single solve() call either:

  - runs a bisection search for μ_max (the largest growth rate the ME-model
    can support given the current constraints), or
  - solves the LP at a fixed μ (when ``mu_fixed`` is provided).

Because each call is a self-contained optimization (not a time-stepped
simulation), the wrapper is a **Step**, not a Process — see the rationale
below.

Why a Step, not a Process
-------------------------
ME-models are constraint-based. There is no internal time integrator — each
call computes a single optimal solution given the current bounds. The
"interval" semantics of a process-bigraph Process don't apply.

A surrounding bigraph that wants TIME DYNAMICS (e.g. dFBA-style integration
of substrate depletion driven by FBA solutions) wires this Step downstream
of a Process that holds the time integration, calling the Step at each
integrator tick to compute new exchange fluxes. The split is the standard
constraint-based pattern (e.g. spatio-flux's dFBA wraps cobra in the same way).

Cell-side interface contract
----------------------------
Designed to conform with
``references/expert/cell_side_interface_contract.md`` in the v2ecoli
multiscale-bioprocess investigation. Inputs (environment + ROS state) and
outputs (exchange fluxes + biomass-mu signal) are the contract surface that
lets this engine serve as a candidate substitute for v2ecoli baseline under
the same reactor composite.

Real bridge fidelity
--------------------
``update()`` calls oxidizeme / cobrame / qminospy directly. There is no
fallback math. If the upstream stack isn't installed, ``update()`` raises
:class:`OxidizeMEMissingError` with a pointer to install instructions. This
is the explicit "guarded real bridge" pattern from the pbg-expert skill.
"""

from __future__ import annotations

from typing import Any

from process_bigraph import Step


# Upstream package surface we depend on. The bridge lazily imports these
# inside ``update()`` so the wrapper file is valid (and shape-tests pass)
# in environments where the upstream stack is absent.
_UPSTREAM_PACKAGES = ("cobra", "cobrame", "ecolime", "qminospy", "oxidizeme")


class OxidizeMEMissingError(RuntimeError):
    """Raised by ``OxidizeMEStep.update`` when the upstream OxidizeME stack
    is not importable. Carries an installation pointer in its message.

    The upstream stack has these *known environmental constraints* as of
    2026-05-22:

    - Python 2.7 (cobrapy 0.5.11 era). Modern Python 3 ports of cobrame
      exist on the ``devel`` branch but are not on PyPI.
    - qMINOS — a 128-bit-precision NLP solver — is **proprietary** and
      obtained by contacting Prof. M. A. Saunders at Stanford University.
      SoPlex (80-bit) is the alternative.
    - ME-models are inherently ill-scaled; only extended-precision solvers
      converge reliably.

    The README documents the install path (DockerHub images, manual qMINOS
    request, SoPlex install).
    """


# Defaults aligned with Yang et al. 2019 and Beulig 2025 OxidizeME runs.
_DEFAULT_ROS_NMOLL = {
    "h2o2_c": 0.2,    # intracellular hydrogen peroxide, nmol/L
    "o2s_c":  0.05,   # intracellular superoxide, nmol/L
}
_DEFAULT_METAL_NMOLL = {
    "fe2_c": 1e4,
    "mn2_c": 1e3,
    "zn2_c": 1e2,
}


class OxidizeMEStep(Step):
    """Real bridge to the OxidizeME ME-model.

    Each ``update()`` call:

      1. Loads the StressME model on first call (heavy — minutes for cobrame
         build-from-scratch; recommended path is a pickled pre-built model
         supplied via ``me_model_pickle_path``).
      2. Applies exchange bounds from ``environment.external_concentrations``
         (mM → mmol/(gDW·h) uptake bounds; the conversion uses a literature
         exchange-flux mapping documented in ``_apply_exchange_bounds``).
      3. Calls ``StressME.substitute_ros`` / ``substitute_metal`` with the
         supplied ROS / metal concentrations.
      4. Calls ``StressME.solve`` (bisection for μ_max if ``mu_fixed`` is
         None; fixed-μ LP otherwise).
      5. Reads the solution: μ, exchange fluxes (mmol/(gDW·h)), full flux
         distribution, computed proteome allocations (translation flux ×
         protein MW per pathway divided by total translation flux), damaged
         proteome fraction.

    Config knobs are the build-time settings (where to load the model from,
    solver verbosity, default ROS levels). Input ports carry the per-call
    state (environment, ROS, mu_fixed override).
    """

    config_schema = {
        "me_model_pickle_path": {
            "_type": "string",
            "_default": "",
            # Path to a pickled StressME model. Strongly recommended over
            # build-from-scratch (which can take minutes and pulls in the
            # full ecolime data tree).
        },
        "build_from_scratch": {
            "_type": "boolean",
            "_default": False,
            # If true and pickle_path is empty, run ``StressME.make_oxidizeme``
            # from a freshly-built ecolime iJL1678b-ME model. Expensive.
        },
        "use_observed": {
            "_type": "boolean",
            "_default": False,
            # Passes ``observed=True`` to StressME constructor — wraps the
            # NLP in ObserveME first.
        },
        "default_ros_nmolL": {
            "_type": "map[float]",
            "_default": _DEFAULT_ROS_NMOLL,
        },
        "default_metal_nmolL": {
            "_type": "map[float]",
            "_default": _DEFAULT_METAL_NMOLL,
        },
        "exchange_id_map": {
            # External-concentration key (cell-side interface contract identifier,
            # e.g. ``GLC[p]``) → upstream ME-model exchange-reaction ID
            # (e.g. ``EX_glc__D_e``). The default below covers the canonical
            # set; override per-experiment.
            "_type": "map[string]",
            "_default": {
                "GLC[p]":             "EX_glc__D_e",
                "OXYGEN-MOLECULE[p]": "EX_o2_e",
                "CARBON-DIOXIDE[p]":  "EX_co2_e",
                "AMMONIUM[p]":        "EX_nh4_e",
                "ACET[p]":            "EX_ac_e",
            },
        },
        "verbose_solver": {"_type": "boolean", "_default": False},
        # Strict mode: if True (default), update() raises OxidizeMEMissingError
        # when the upstream stack is absent. If False, returns an empty
        # solution dict — useful only for dashboard scaffolding.
        "strict": {"_type": "boolean", "_default": True},
    }

    def __init__(self, config: dict | None = None, core: Any = None) -> None:
        super().__init__(config=config, core=core)
        self._stress_me: Any | None = None  # StressME instance (lazy)
        self._build_attempted = False
        self._build_error: Exception | None = None
        self._upstream_missing: list[str] = []  # populated by _import_upstream

    # ─── Ports ───────────────────────────────────────────────────────────

    def inputs(self) -> dict[str, str]:
        # All inputs are state a sibling Process can sensibly write to.
        # See cell_side_interface_contract.md.
        return {
            "external_concentrations": "map[float]",   # mM, keyed by EcoCyc-style IDs
            "ros_concentrations":      "map[float]",   # nmol/L, keys: h2o2_c, o2s_c, ...
            "metal_concentrations":    "map[float]",   # nmol/L, keys: fe2_c, mn2_c, zn2_c
            "mu_fixed":                "maybe[float]", # 1/h; None ⇒ bisection
            "atpm_mmol_per_gDW_per_h": "maybe[float]", # ATP maintenance reaction lower bound override
        }

    def outputs(self) -> dict[str, str]:
        # Solution snapshot. ME-model "results" are absolute (not deltas),
        # so we use ``overwrite[float]`` for scalar setpoints (mu, damaged
        # fraction, status) per the cell-side interface contract's
        # "replace semantics for tool-internal absolute readings" guidance.
        return {
            "mu":                          "overwrite[float]",   # optimal specific growth rate, 1/h
            "external_exchange_fluxes":    "map[float]",         # mmol/(gDW·h), keyed by EcoCyc-style IDs
            "proteome_allocations":        "map[float]",         # fraction per pathway/sector
            "damaged_proteome_fraction":   "overwrite[float]",   # 0..1
            "solver_status":               "overwrite[string]",  # 'optimal' | 'infeasible' | 'unbounded' | ...
            "wall_time_s":                 "overwrite[float]",   # solver wall time per call
        }

    def initial_state(self) -> dict[str, Any]:
        return {
            "external_concentrations": {},
            "ros_concentrations":      dict(self.config["default_ros_nmolL"]),
            "metal_concentrations":    dict(self.config["default_metal_nmolL"]),
            "mu_fixed":                None,
            "atpm_mmol_per_gDW_per_h": None,
        }

    # ─── Update ──────────────────────────────────────────────────────────

    def update(self, state: dict[str, Any]) -> dict[str, Any]:
        # Lazy import the upstream stack. The whole real-bridge fidelity
        # lives below this line.
        upstream = self._import_upstream()
        if upstream is None:
            if self.config["strict"]:
                raise OxidizeMEMissingError(self._missing_message())
            # Non-strict scaffolding fallback: return a labelled empty solution
            # so the dashboard can wire the Step without crashing. Numbers are
            # NaN where we can't fake science honestly.
            import math
            return {
                "mu":                        float("nan"),
                "external_exchange_fluxes":  {},
                "proteome_allocations":      {},
                "damaged_proteome_fraction": float("nan"),
                "solver_status":             "upstream_missing",
                "wall_time_s":               0.0,
            }

        if self._stress_me is None:
            self._build_stress_me(upstream)

        return self._run_one_solve(state, upstream)

    # ─── Real-bridge internals ──────────────────────────────────────────

    def _import_upstream(self) -> dict[str, Any] | None:
        """Lazily import the upstream stack. Returns a dict of module handles
        if all packages are present; otherwise records the missing names and
        returns None. Idempotent across calls.
        """
        if self._upstream_missing:
            return None
        modules: dict[str, Any] = {}
        for name in _UPSTREAM_PACKAGES:
            try:
                modules[name] = __import__(name)
            except ImportError as exc:
                self._upstream_missing.append(f"{name} ({exc})")
        if self._upstream_missing:
            return None
        # Additionally pull StressME by attribute access (oxidizeme exposes it
        # via oxidizeme.model.StressME).
        from oxidizeme.model import StressME       # noqa: WPS433  (lazy)
        from qminospy.me1 import ME_NLP1           # noqa: WPS433
        modules["StressME"] = StressME
        modules["ME_NLP1"]  = ME_NLP1
        return modules

    def _missing_message(self) -> str:
        missing = ", ".join(self._upstream_missing) or "unknown"
        return (
            f"OxidizeMEStep cannot run: upstream packages are not importable "
            f"({missing}). The upstream stack is Python 2.7 / cobrapy 0.5.11 "
            f"era; qMINOS is proprietary (request from Prof. M. A. Saunders, "
            f"Stanford). Install paths:\n"
            f"  1. Recommended — pull SBRG's Docker image with the full stack:\n"
            f"     https://hub.docker.com/r/sbrg/cobrame\n"
            f"  2. Source builds — see README.md ('Real bridge — install obstacle').\n"
            f"  3. SoPlex (80-bit, free for academic use) as a qMINOS alternative.\n"
            f"To wire OxidizeMEStep into a dashboard composite WITHOUT running the "
            f"solver, set ``config.strict = False`` (returns a labelled empty solution; "
            f"NEVER assert numerical fidelity on those outputs)."
        )

    def _build_stress_me(self, upstream: dict[str, Any]) -> None:
        """Load (or build) the underlying StressME model. Heavy operation;
        memoised on the Step instance.
        """
        if self._build_attempted:
            if self._build_error is not None:
                raise self._build_error
            return
        self._build_attempted = True
        try:
            pkl = self.config["me_model_pickle_path"]
            if pkl:
                import pickle
                with open(pkl, "rb") as fh:
                    me_nlp = pickle.load(fh)
                stress_me = upstream["StressME"](
                    me_nlp,
                    observed=self.config["use_observed"],
                )
            elif self.config["build_from_scratch"]:
                stress_me = self._build_from_scratch(upstream)
            else:
                raise OxidizeMEMissingError(
                    "OxidizeMEStep needs either ``me_model_pickle_path`` or "
                    "``build_from_scratch=true`` to load the StressME model. "
                    "Build-from-scratch is expensive (~minutes); pickling a "
                    "pre-built model is the recommended path."
                )
            self._stress_me = stress_me
        except Exception as exc:
            self._build_error = exc
            raise

    def _build_from_scratch(self, upstream: dict[str, Any]) -> Any:
        """Build the StressME model from ecolime + apply ``make_oxidizeme``.
        Heavy. See README "Real bridge — building from scratch".
        """
        from ecolime import build_ME_model
        # ecolime.build_ME_model.return_me_model() is the canonical entry
        # point in the 0.0.9 era; older versions named it differently.
        me_model = build_ME_model.return_me_model()
        me_nlp = upstream["ME_NLP1"](me_model, growth_key="mu")
        stress_me = upstream["StressME"](
            me_nlp,
            observed=self.config["use_observed"],
        )
        stress_me.make_oxidizeme(
            force_damage=True,
            extra_dilution=True,
        )
        return stress_me

    def _run_one_solve(
        self,
        state: dict[str, Any],
        upstream: dict[str, Any],
    ) -> dict[str, Any]:
        """Apply the per-call constraints, solve, harvest the solution."""
        import time
        stress_me = self._stress_me

        # 1) Exchange bounds from external_concentrations (mM → uptake bound)
        self._apply_exchange_bounds(stress_me, state.get("external_concentrations", {}))

        # 2) ATPM override if supplied
        atpm = state.get("atpm_mmol_per_gDW_per_h")
        if atpm is not None and atpm > 0:
            # Standard cobrapy lower-bound override for the ATP maintenance
            # reaction (cobrame variants: 'ATPM' or 'ATPM_FWD').
            for rid in ("ATPM", "ATPM_FWD"):
                try:
                    stress_me.me.reactions.get_by_id(rid).lower_bound = float(atpm)
                    break
                except (KeyError, AttributeError):
                    continue

        # 3) ROS + metal substitutions
        ros = state.get("ros_concentrations") or {}
        if ros:
            stress_me.substitute_ros(stress_me.me_nlp, subs_dict=ros)
        metals = state.get("metal_concentrations") or {}
        if metals:
            stress_me.substitute_metal(stress_me.me_nlp, subs_dict=metals)

        # 4) Solve
        t0 = time.perf_counter()
        mu_fixed = state.get("mu_fixed")
        try:
            mu_opt = stress_me.solve(mu_fixed=mu_fixed)
            status = "optimal"
        except Exception as exc:
            status = f"error:{type(exc).__name__}"
            mu_opt = float("nan")
        wall = time.perf_counter() - t0

        # 5) Harvest solution
        return {
            "mu":                        float(mu_opt) if mu_opt == mu_opt else float("nan"),
            "external_exchange_fluxes":  self._harvest_exchange_fluxes(stress_me),
            "proteome_allocations":      self._compute_proteome_allocations(stress_me),
            "damaged_proteome_fraction": self._compute_damaged_proteome(stress_me),
            "solver_status":             status,
            "wall_time_s":               wall,
        }

    def _apply_exchange_bounds(
        self,
        stress_me: Any,
        external_concentrations_mM: dict[str, float],
    ) -> None:
        """Map external_concentrations (mM) into exchange-reaction lower bounds.

        Convention: a higher external concentration relaxes the uptake bound
        (more negative lower bound on the exchange reaction). The mapping
        below is a textbook Michaelis-style soft cap; the spec PR (in the
        v2ecoli mbp investigation) replaces this with a per-substrate
        calibrated form derived from the cell-side interface contract.
        """
        # K-half and v-max defaults; replaceable per-substrate via config in
        # a future iteration. These are placeholder calibration values
        # ANCHORED IN literature ranges for E. coli (μmol/(gDW·h)·OD-scale).
        K_HALF_MM = 0.1
        V_MAX_MMOL_GDW_H = {
            "EX_glc__D_e":  10.0,
            "EX_o2_e":      20.0,
            "EX_nh4_e":     10.0,
            "EX_co2_e":     50.0,  # excretion, upper bound
            "EX_ac_e":      20.0,
        }
        exchange_id_map = self.config["exchange_id_map"]
        for ext_key, conc_mM in external_concentrations_mM.items():
            rxn_id = exchange_id_map.get(ext_key)
            if not rxn_id:
                continue
            try:
                rxn = stress_me.me.reactions.get_by_id(rxn_id)
            except (KeyError, AttributeError):
                continue
            v_max = V_MAX_MMOL_GDW_H.get(rxn_id, 10.0)
            # Saturating uptake limit (Monod-style)
            v = -v_max * (conc_mM / (K_HALF_MM + conc_mM))
            rxn.lower_bound = float(v)

    def _harvest_exchange_fluxes(self, stress_me: Any) -> dict[str, float]:
        """Read exchange-flux solution back from stress_me.me.solution."""
        exchange_id_map = self.config["exchange_id_map"]
        out: dict[str, float] = {}
        try:
            x = stress_me.me.solution.x_dict
        except (AttributeError, TypeError):
            return out
        # Invert exchange_id_map so we can emit on the cell-side keys
        rxn_to_ext = {v: k for k, v in exchange_id_map.items()}
        for rxn_id, flux in x.items():
            ext_key = rxn_to_ext.get(rxn_id)
            if ext_key is not None:
                out[ext_key] = float(flux)
        return out

    def _compute_proteome_allocations(self, stress_me: Any) -> dict[str, float]:
        """Replicate Yang 2019 / Beulig 2025 proteome-allocation formula:

            % proteome_i = (Σ_{r∈pathway_i} mw_r × v_translation_r)
                         / (Σ_{r∈all}        mw_r × v_translation_r)

        Returns fractions per named pathway (glycolysis, PPP, TCA cycle, and
        a catch-all 'other'). Per Beulig 2025's Methods.
        """
        out: dict[str, float] = {}
        try:
            x = stress_me.me.solution.x_dict
            translation_rxns = [r for r in stress_me.me.reactions if r.id.startswith("translation_")]
        except (AttributeError, TypeError):
            return out

        # Map locus → molecular weight from the underlying ME-model.
        # cobrame's TranslationReaction carries an associated Protein object
        # which has a molecular_weight attribute (Da).
        total = 0.0
        per_pathway: dict[str, float] = {"glycolysis": 0.0, "ppp": 0.0, "tca": 0.0, "other": 0.0}
        # Pathway membership (rough; expand via config in a future iteration)
        PATHWAY_KEYWORDS = {
            "glycolysis": ("PGI", "PFK", "FBA", "TPI", "GAPD", "PGK", "PGM", "ENO", "PYK"),
            "ppp":        ("G6PDH2r", "PGL", "GND", "RPI", "RPE", "TKT", "TALA"),
            "tca":        ("CS", "ACONTa", "ACONTb", "ICDHyr", "AKGDH", "SUCOAS", "SUCDi", "FUM", "MDH"),
        }
        # Build locus → pathway mapping by inspecting reaction associations
        # (best-effort; cobrame's translation rxn IDs include the locus tag).
        for trans_rxn in translation_rxns:
            v = float(x.get(trans_rxn.id, 0.0))
            try:
                mw = trans_rxn.protein.molecular_weight  # Da
            except AttributeError:
                continue
            contribution = v * mw
            total += contribution
            # Classify by linked metabolic reactions
            classified = "other"
            try:
                gene = trans_rxn.id.replace("translation_", "")
                # Look up reactions catalysed by this gene's protein
                related_rxn_ids = {
                    rxn.id for rxn in stress_me.me.reactions
                    if hasattr(rxn, "complex_data")
                    and rxn.complex_data is not None
                    and gene in str(rxn.complex_data)
                }
                for pw, keywords in PATHWAY_KEYWORDS.items():
                    if any(kw in rid for rid in related_rxn_ids for kw in keywords):
                        classified = pw
                        break
            except Exception:
                pass
            per_pathway[classified] += contribution

        if total > 0:
            out = {pw: float(per_pathway[pw] / total) for pw in per_pathway}
        return out

    def _compute_damaged_proteome(self, stress_me: Any) -> float:
        """Replicate Yang 2019 / Beulig 2025 damaged-proteome formula:

            % damaged_proteome = (Σ_{r∈damage_*} mw_r × v_ComplexFormation_r)
                               / (Σ_{r∈translation} mw_r × v_translation_r)
        """
        try:
            x = stress_me.me.solution.x_dict
        except (AttributeError, TypeError):
            return float("nan")

        translation_total = 0.0
        for r in stress_me.me.reactions:
            if r.id.startswith("translation_"):
                try:
                    mw = r.protein.molecular_weight
                    translation_total += float(x.get(r.id, 0.0)) * mw
                except AttributeError:
                    continue

        damage_total = 0.0
        for r in stress_me.me.reactions:
            if r.id.startswith("damage_") and r.id.endswith("ComplexFormation"):
                try:
                    mw = r.protein.molecular_weight
                    damage_total += float(x.get(r.id, 0.0)) * mw
                except AttributeError:
                    continue

        if translation_total <= 0:
            return float("nan")
        return float(damage_total / translation_total)
