"""OxidizeME steady-state composite.

The ME-model is constraint-based — there is no time integration internal to
OxidizeME. The natural "composite" wired around it is a single-shot
steady-state optimization driven by an environment + ROS context, with
an emitter capturing the solution.

For a time-dynamic context (dFBA-style — substrate depletion → re-optimize
periodically), wire OxidizeMEStep DOWNSTREAM of a Process that holds the
time integration; that Process drives external_concentrations / ros each
tick and reads back μ + exchange_fluxes. See the v2ecoli multiscale-bioprocess
investigation (mbp-03-bird-reactor-coupling) for a worked example.
"""

from viva_superpowers.composite_generator import composite_generator


@composite_generator(
    name="oxidizeme_steady_state",
    description="Single-shot OxidizeME ME-model solve at given environment + ROS context.",
    parameters={
        "me_model_pickle_path": {
            "type": "string",
            "default": "",
            "description": "Path to pickled StressME model (strongly recommended)",
        },
        "glucose_mM": {
            "type": "float",
            "default": 22.2,    # ~4 g/L M9-glucose
            "description": "External glucose concentration (mM)",
        },
        "oxygen_mM": {
            "type": "float",
            "default": 0.21,    # near-saturation aerobic
            "description": "External dissolved O2 (mM)",
        },
        "ammonium_mM": {
            "type": "float",
            "default": 19.0,    # M9 ~1 g/L NH4Cl
            "description": "External ammonium (mM)",
        },
        "h2o2_c_nmolL": {
            "type": "float",
            "default": 0.2,
            "description": "Intracellular H2O2 (nmol/L)",
        },
        "o2s_c_nmolL": {
            "type": "float",
            "default": 0.05,
            "description": "Intracellular superoxide (nmol/L)",
        },
        "mu_fixed": {
            "type": "float",
            "default": -1.0,
            "description": "Fixed μ (1/h); -1 means bisect for μ_max",
        },
        "strict": {
            "type": "boolean",
            "default": False,
            "description": "If false, returns NaN solution when OxidizeME stack missing (dashboard-friendly)",
        },
    },
)
def oxidizeme_steady_state(
    core=None,
    *,
    me_model_pickle_path: str = "",
    glucose_mM: float = 22.2,
    oxygen_mM: float = 0.21,
    ammonium_mM: float = 19.0,
    h2o2_c_nmolL: float = 0.2,
    o2s_c_nmolL: float = 0.05,
    mu_fixed: float = -1.0,
    strict: bool = False,
):
    """Build the steady-state composite document."""
    mu = None if mu_fixed < 0 else mu_fixed
    return {
        "oxidizeme": {
            "_type": "step",
            "address": "local:OxidizeMEStep",
            "config": {
                "me_model_pickle_path": me_model_pickle_path,
                "strict": strict,
            },
            "inputs": {
                "external_concentrations": ["stores", "external_concentrations"],
                "ros_concentrations":      ["stores", "ros_concentrations"],
                "metal_concentrations":    ["stores", "metal_concentrations"],
                "mu_fixed":                ["stores", "mu_fixed"],
                "atpm_mmol_per_gDW_per_h": ["stores", "atpm"],
            },
            "outputs": {
                "mu":                          ["stores", "mu"],
                "external_exchange_fluxes":    ["stores", "exchange_fluxes"],
                "proteome_allocations":        ["stores", "proteome_allocations"],
                "damaged_proteome_fraction":   ["stores", "damaged_proteome"],
                "solver_status":               ["stores", "solver_status"],
                "wall_time_s":                 ["stores", "wall_time_s"],
            },
        },
        "stores": {
            "external_concentrations": {
                "GLC[p]":             float(glucose_mM),
                "OXYGEN-MOLECULE[p]": float(oxygen_mM),
                "AMMONIUM[p]":        float(ammonium_mM),
            },
            "ros_concentrations": {
                "h2o2_c": float(h2o2_c_nmolL),
                "o2s_c":  float(o2s_c_nmolL),
            },
            "metal_concentrations": {
                "fe2_c": 1e4,
                "mn2_c": 1e3,
                "zn2_c": 1e2,
            },
            "mu_fixed": mu,
            "atpm":     None,
            "mu":                 0.0,
            "exchange_fluxes":    {},
            "proteome_allocations": {},
            "damaged_proteome":   0.0,
            "solver_status":      "",
            "wall_time_s":        0.0,
        },
        "emitter": {
            "_type": "step",
            "address": "local:RAMEmitter",
            "config": {
                "emit": {
                    "mu":                 "float",
                    "exchange_fluxes":    "map[float]",
                    "proteome_allocations": "map[float]",
                    "damaged_proteome":   "float",
                    "solver_status":      "string",
                    "wall_time_s":        "float",
                },
            },
            "inputs": {
                "mu":                 ["stores", "mu"],
                "exchange_fluxes":    ["stores", "exchange_fluxes"],
                "proteome_allocations": ["stores", "proteome_allocations"],
                "damaged_proteome":   ["stores", "damaged_proteome"],
                "solver_status":      ["stores", "solver_status"],
                "wall_time_s":        ["stores", "wall_time_s"],
            },
        },
    }
