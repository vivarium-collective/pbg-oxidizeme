"""Demo report for pbg-oxidizeme.

This demo is honest about what's running and what isn't.

WHAT RUNS:
- iML1515 (E. coli K-12 MG1655 genome-scale M-model from BiGG, the substrate
  OxidizeME extends with macromolecular expression + ROS damage) — solved at
  three operating conditions via cobrapy, on Python 3, no proprietary deps.
- Real Beulig 2025 (mSystems) fermentation trajectories, loaded from the
  digitized CSVs at v2ecoli/.../palsson-2025-supp.

WHAT DOESN'T RUN HERE:
- OxidizeME itself (the ME-model layer with expression + ROS) — the upstream
  stack is Python 2.7 / cobrapy 0.5.11 / qMINOS (proprietary). The wrapper
  raises OxidizeMEMissingError in strict mode and returns labelled empty
  solutions in non-strict mode. See README "Real bridge — install obstacle".

WHAT THE DEMO SHOWS:
- The wrapper's port surface + composite architecture.
- iML1515 as a M-model preview: flux distributions, μ vs O2 / glucose, the
  acetate-overflow signature OxidizeME inherits.
- Beulig 2025 measured trajectories: OD, glucose, growth rate, OTR/CTR,
  byproducts across the batch and fed-batch phases of real reactors.
- Side-by-side: M-model predicted μ vs Beulig measured μ at matched conditions.
- A "What OxidizeME adds on top of iML1515" panel making the upgrade concrete.

Output: demo/report.html
"""

from __future__ import annotations

import base64
import json
import math
import os
import sys
import time
import webbrowser
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = REPO_ROOT / "demo"
REPORT_PATH = OUT_DIR / "report.html"

# Beulig 2025 supplementary data — already mirrored into the v2ecoli mbp
# investigation. Use a configurable path so this demo runs even when the
# v2ecoli checkout isn't where we expect.
DEFAULT_BEULIG_DIR = REPO_ROOT / "references" / "papers" / "palsson-2025-supp"
BEULIG_DIR = Path(os.environ.get("BEULIG_SUPP", str(DEFAULT_BEULIG_DIR)))


# ─── iML1515 (the M-model OxidizeME extends) ────────────────────────────────


def _disable_o2_uptake(model, bound: float = -1000.0) -> None:
    """Override iML1515's EX_o2_e lower bound to limit O2 availability.
    Less-negative bound = tighter O2 supply."""
    rxn = model.reactions.get_by_id("EX_o2_e")
    rxn.lower_bound = bound


def _disable_glc_uptake(model, bound: float = -10.0) -> None:
    """Override iML1515's EX_glc__D_e lower bound to limit glucose availability."""
    rxn = model.reactions.get_by_id("EX_glc__D_e")
    rxn.lower_bound = bound


def run_iml1515_configs() -> list[dict[str, Any]]:
    """Run iML1515 at three conditions and return structured results.

    Each result records the config, μ, key exchange fluxes (glucose, O2, CO2,
    acetate, ethanol, formate, lactate, succinate), and wall time. These are
    REAL FBA solutions, not placeholders.
    """
    import cobra
    print("[iml1515] loading model…")
    model = cobra.io.load_model("iML1515")

    # Sanity: cache the originals so each config starts from the default bounds
    orig_o2  = model.reactions.get_by_id("EX_o2_e").lower_bound
    orig_glc = model.reactions.get_by_id("EX_glc__D_e").lower_bound

    configs = [
        {
            "id": "aerobic-baseline",
            "title": "Aerobic baseline (M9 + glucose, ample O2)",
            "subtitle": "iML1515 with default exchange bounds — the reference μ_max",
            "description": (
                "iML1515's default M9-glucose-aerobic configuration. Establishes "
                "the M-model reference μ; OxidizeME's ME-layer will modulate this "
                "via expression costs and ROS damage in the full bridge."
            ),
            "o2_bound":  -1000.0,
            "glc_bound": -10.0,
            "accent": "#10b981",
        },
        {
            "id": "low-o2",
            "title": "O2-limited (mimics fed-batch onset)",
            "subtitle": "Tight EX_o2_e bound — overflow regime emerges",
            "description": (
                "Tightening the O2 uptake bound to −5 mmol/(gDW·h) drops μ and pushes "
                "carbon into fermentation byproducts (overflow). This is the regime "
                "Beulig 2025 sees at high cell density as kLa stops keeping up with "
                "the population's O2 demand."
            ),
            "o2_bound":  -5.0,
            "glc_bound": -10.0,
            "accent": "#f59e0b",
        },
        {
            "id": "glucose-limited",
            "title": "Glucose-limited",
            "subtitle": "Tight EX_glc__D_e bound — feed-limited fed-batch analogue",
            "description": (
                "Tightening the glucose uptake bound to −5 mmol/(gDW·h) — exponential-"
                "feed conditions where glucose arrives slower than maximal demand. μ "
                "drops in lockstep with the bound; byproduct excretion falls."
            ),
            "o2_bound":  -1000.0,
            "glc_bound": -5.0,
            "accent": "#3b82f6",
        },
    ]

    exchange_keys = [
        "EX_glc__D_e", "EX_o2_e", "EX_co2_e", "EX_ac_e", "EX_etoh_e",
        "EX_for_e", "EX_lac__D_e", "EX_succ_e", "EX_nh4_e", "BIOMASS_Ec_iML1515_core_75p37M",
    ]

    results: list[dict[str, Any]] = []
    for cfg in configs:
        # Reset to defaults each time
        model.reactions.get_by_id("EX_o2_e").lower_bound  = orig_o2
        model.reactions.get_by_id("EX_glc__D_e").lower_bound = orig_glc
        _disable_o2_uptake(model, bound=cfg["o2_bound"])
        _disable_glc_uptake(model, bound=cfg["glc_bound"])

        t0 = time.perf_counter()
        sol = model.optimize()
        wall = time.perf_counter() - t0

        # Capture exchange fluxes via the cached solution.fluxes Series
        fluxes = {}
        for rid in exchange_keys:
            try:
                fluxes[rid] = float(sol.fluxes[rid])
            except (KeyError, AttributeError):
                fluxes[rid] = float("nan")
        # Top-10 absolute flux reactions for the long table
        try:
            top_fluxes = sol.fluxes.abs().nlargest(15).to_dict()
            top_fluxes = {k: float(sol.fluxes[k]) for k in top_fluxes}
        except Exception:
            top_fluxes = {}

        results.append({
            "config": cfg,
            "mu": float(sol.objective_value) if sol.objective_value is not None else float("nan"),
            "status": str(sol.status),
            "wall_s": wall,
            "exchanges": fluxes,
            "top_fluxes": top_fluxes,
        })
    return results


# ─── Beulig 2025 measured trajectories ──────────────────────────────────────


def _read_csv_dict(path: Path) -> list[dict[str, Any]]:
    import csv
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def _to_float(s: Any) -> float | None:
    try:
        v = float(s)
        if math.isnan(v):
            return None
        return v
    except (TypeError, ValueError):
        return None


def load_beulig_trajectories(reactor_ids: list[str] | None = None) -> dict[str, Any]:
    """Read process_summary_interpol.csv and return per-reactor time-series.

    Returns a dict {
      'reactor_ids': [...],
      'series': {reactor_id: {'time_h': [...], 'OD': [...], 'glucose_mmolL': [...], ...}, ...},
      'available': bool — False if BEULIG_DIR missing
    }
    """
    summary = BEULIG_DIR / "process_summary_interpol.csv"
    if not summary.exists():
        return {"reactor_ids": [], "series": {}, "available": False, "missing_path": str(summary)}

    rows = _read_csv_dict(summary)
    # Reactors with the most observations are best for charts
    counts: dict[str, int] = {}
    for r in rows:
        rid = r.get("reactor_id", "").strip()
        if rid:
            counts[rid] = counts.get(rid, 0) + 1

    # Pick the 6 reactors with the densest sampling unless caller specifies.
    if reactor_ids is None:
        top = sorted(counts.items(), key=lambda kv: -kv[1])[:6]
        reactor_ids = [rid for rid, _ in top]

    series: dict[str, dict[str, list[float]]] = {rid: {
        "time_h": [], "OD": [], "glucose_mmolL": [], "growth_rate_per_h": [],
        "OTR_molh": [], "CTR_molh": [], "feed_rate_Lph": [],
        "ethanol_mmolL": [], "formate_mmolL": [],
    } for rid in reactor_ids}

    for r in rows:
        rid = r.get("reactor_id", "").strip()
        if rid not in series:
            continue
        t = _to_float(r.get("fed-batch time [h]"))
        if t is None:
            continue
        series[rid]["time_h"].append(t)
        series[rid]["OD"].append(_to_float(r.get("OD [-]")) or float("nan"))
        series[rid]["glucose_mmolL"].append(_to_float(r.get("D-glucose [mmol/L]")) or float("nan"))
        series[rid]["growth_rate_per_h"].append(_to_float(r.get("growth rate [1/h]")) or float("nan"))
        series[rid]["OTR_molh"].append(_to_float(r.get("OTR [mol/h]")) or float("nan"))
        series[rid]["CTR_molh"].append(_to_float(r.get("CTR [mol/h]")) or float("nan"))
        series[rid]["feed_rate_Lph"].append(_to_float(r.get("feed rate [L/h]")) or float("nan"))
        series[rid]["ethanol_mmolL"].append(_to_float(r.get("ethanol [mmol/L]")) or float("nan"))
        series[rid]["formate_mmolL"].append(_to_float(r.get("formate [mmol/L]")) or float("nan"))

    # Sort each reactor's points by time
    for rid in reactor_ids:
        idx = sorted(range(len(series[rid]["time_h"])), key=lambda i: series[rid]["time_h"][i])
        for k, v in series[rid].items():
            series[rid][k] = [v[i] for i in idx]

    return {"reactor_ids": reactor_ids, "series": series, "available": True}


# ─── Bigraph PNG ─────────────────────────────────────────────────────────────


def render_bigraph_png() -> str | None:
    try:
        from bigraph_viz import plot_bigraph
    except ImportError:
        return None
    from pbg_oxidizeme.composites import oxidizeme_steady_state

    doc = oxidizeme_steady_state(core=None, strict=False)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    try:
        plot_bigraph(
            state=doc,
            out_dir=str(OUT_DIR),
            filename="bigraph",
            file_format="png",
            remove_process_place_edges=True,
            rankdir="LR",
            port_labels=False,
            dpi="150",
        )
    except Exception as exc:
        print(f"bigraph-viz render failed: {exc}", file=sys.stderr)
        return None
    png = OUT_DIR / "bigraph.png"
    if not png.exists():
        return None
    with open(png, "rb") as fh:
        return "data:image/png;base64," + base64.b64encode(fh.read()).decode()


# ─── HTML rendering ──────────────────────────────────────────────────────────


def _fmt_flux(v: float) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "NaN"
    if abs(v) < 1e-6:
        return f"{v:.2e}"
    return f"{v:.4g}"


def _flux_row(rid: str, v: float) -> str:
    direction = "→" if v >= 0 else "←"
    cls = "flux-up" if v >= 0 else "flux-down"
    return f'<tr><td><code>{rid}</code></td><td class="{cls}">{direction} {_fmt_flux(abs(v))}</td></tr>'


def render_iml1515_card(result: dict[str, Any]) -> str:
    cfg = result["config"]
    flux_rows = "".join(_flux_row(rid, v) for rid, v in result["exchanges"].items())
    top_rows = "".join(_flux_row(rid, v) for rid, v in result["top_fluxes"].items())
    return f"""
<section class="card" style="border-left-color:{cfg['accent']}">
  <header>
    <h3>{cfg['title']}</h3>
    <div class="subtitle">{cfg['subtitle']}</div>
    <div class="badges">
      <span class="badge ok">μ = {result['mu']:.4f} 1/h</span>
      <span class="badge">status: {result['status']}</span>
      <span class="badge muted">EX_o2_e ≥ {cfg['o2_bound']:g}, EX_glc__D_e ≥ {cfg['glc_bound']:g}</span>
      <span class="badge muted">{result['wall_s']*1000:.0f} ms</span>
    </div>
  </header>
  <p>{cfg['description']}</p>
  <div class="two-col">
    <div>
      <h4>Key exchange fluxes (mmol/(gDW·h))</h4>
      <table>{flux_rows}</table>
    </div>
    <div>
      <h4>Top |fluxes| across the network</h4>
      <table>{top_rows}</table>
    </div>
  </div>
</section>
"""


def render_beulig_section(beulig: dict[str, Any]) -> str:
    if not beulig["available"]:
        return f"""
<section class="card warn">
  <h3>Beulig 2025 trajectories</h3>
  <p>Data not available — set <code>BEULIG_SUPP</code> env var to point at
  <code>references/papers/palsson-2025-supp/</code>. Looked for:
  <code>{beulig.get('missing_path','?')}</code></p>
</section>"""

    # Embed the per-reactor series as JSON so Plotly can draw client-side
    payload = json.dumps(beulig["series"], default=lambda o: float("nan") if isinstance(o, float) and math.isnan(o) else o)
    n = len(beulig["reactor_ids"])
    return f"""
<section class="card" style="border-left-color:#7c3aed">
  <header>
    <h3>Beulig 2025 measured fermentation trajectories</h3>
    <div class="subtitle">Real data from {n} reactors in <code>process_summary_interpol.csv</code> (mirrored from <a href="https://github.com/febedtu/hd_ecoli" target="_blank">febedtu/hd_ecoli</a>, MIT)</div>
  </header>
  <p>Each reactor has time-aligned OD, glucose, growth rate, OTR/CTR, byproducts.
  Negative fed-batch time = batch phase (the slice mbp-05 scopes to until pbg-bioreactordesign
  gains fed-batch). Positive = fed-batch phase.</p>

  <div class="grid-2">
    <div id="chart-od"      class="plotly"></div>
    <div id="chart-glucose" class="plotly"></div>
    <div id="chart-mu"      class="plotly"></div>
    <div id="chart-otr"     class="plotly"></div>
    <div id="chart-feed"    class="plotly"></div>
    <div id="chart-ethanol" class="plotly"></div>
  </div>

<script>
const BEULIG = {payload};
const REACTORS = Object.keys(BEULIG);
function trace(rid, ykey) {{
  const s = BEULIG[rid];
  return {{x: s.time_h, y: s[ykey], mode: "lines+markers", name: rid,
           marker: {{size: 4}}, line: {{width: 1.5}}}};
}}
function make(divId, ykey, ytitle) {{
  const data = REACTORS.map(r => trace(r, ykey));
  Plotly.newPlot(divId, data, {{
    title: ytitle, margin: {{l: 60, r: 10, t: 40, b: 50}}, height: 280,
    xaxis: {{title: "fed-batch time [h]"}}, yaxis: {{title: ytitle}},
    legend: {{font: {{size: 9}}, orientation: "h", y: -0.25}},
  }}, {{displaylogo: false, displayModeBar: false}});
}}
make("chart-od",      "OD",                 "OD600");
make("chart-glucose", "glucose_mmolL",      "D-glucose [mmol/L]");
make("chart-mu",      "growth_rate_per_h",  "growth rate [1/h]");
make("chart-otr",     "OTR_molh",           "OTR [mol/h]");
make("chart-feed",    "feed_rate_Lph",      "feed rate [L/h]");
make("chart-ethanol", "ethanol_mmolL",      "ethanol [mmol/L]");
</script>
</section>
"""


def render_iml1515_vs_beulig(iml: list[dict[str, Any]], beulig: dict[str, Any]) -> str:
    """Side-by-side comparison panel — iML1515-predicted μ at the demo's
    three O2 / glucose bounds, plotted against the range of growth rates
    Beulig actually measured.
    """
    iml_mu = [r["mu"] for r in iml]
    iml_labels = [r["config"]["id"] for r in iml]

    # Collect Beulig measured μ across all reactors and times for context.
    measured_mu: list[float] = []
    if beulig["available"]:
        for s in beulig["series"].values():
            for v in s["growth_rate_per_h"]:
                if isinstance(v, (int, float)) and not math.isnan(v):
                    measured_mu.append(float(v))

    payload = json.dumps({"iml_mu": iml_mu, "iml_labels": iml_labels, "measured": measured_mu})
    return f"""
<section class="card" style="border-left-color:#0ea5e9">
  <header>
    <h3>iML1515-predicted μ vs Beulig 2025 measured μ</h3>
    <div class="subtitle">M-model bounds-driven μ versus the empirical distribution from {len(measured_mu)} measurements</div>
  </header>
  <p>iML1515 at the three demo bounds gives the bars below. The histogram is
  growth rates Beulig actually measured across all reactors and time points
  in <code>process_summary_interpol.csv</code>. OxidizeME adds the
  expression-cost / damage layer that modulates this M-model μ.</p>
  <div id="chart-compare" class="plotly" style="height:340px"></div>
<script>
const CMP = {payload};
Plotly.newPlot("chart-compare", [
  {{
    x: CMP.measured, type: "histogram", name: "Beulig measured μ",
    marker: {{color: "rgba(124,58,237,0.55)"}}, opacity: 0.7,
    histnorm: "probability density", xbins: {{size: 0.02}},
  }},
  {{
    x: CMP.iml_mu, y: [0.4]*CMP.iml_mu.length, mode: "markers+text",
    text: CMP.iml_labels, textposition: "top center", name: "iML1515 demo configs",
    marker: {{size: 16, color: "#0ea5e9", line: {{color:"#0369a1", width:2}}}}, yaxis: "y2",
  }},
], {{
  margin: {{l: 60, r: 60, t: 30, b: 50}},
  xaxis: {{title: "specific growth rate μ [1/h]", range: [0, 1.2]}},
  yaxis: {{title: "Beulig density"}},
  yaxis2: {{overlaying: "y", side: "right", showticklabels: false, range: [0, 1]}},
  legend: {{orientation: "h", y: 1.1}},
}}, {{displaylogo: false, displayModeBar: false}});
</script>
</section>
"""


def render_what_oxidizeme_adds() -> str:
    return """
<section class="card" style="border-left-color:#ec4899">
  <header><h3>What OxidizeME adds on top of iML1515 (the bridge target)</h3></header>
  <p>iML1515 is an M-model — pure flux balance over the metabolic network. OxidizeME
  (Yang et al. 2019 PNAS) replaces the metabolic-only constraint set with a
  <strong>metabolism + macromolecular expression (ME)</strong> model that also tracks:</p>
  <ul>
    <li><strong>Protein synthesis as explicit reactions</strong> — every enzyme has a
        translation reaction; expression flux ≈ enzyme abundance × dilution rate;
        the ME-LP must allocate enough proteome to carry the metabolic flux.</li>
    <li><strong>Ribosome allocation</strong> — translation flux is bounded by ribosome
        availability, which itself is a translation product.</li>
    <li><strong>ROS damage and repair reactions</strong> — H₂O₂ / O₂⁻ damage Fe-S clusters
        and metal-dependent enzymes; repair complexes (IscU, alt-metallation, etc.)
        compete for proteome. Damaged enzymes are diluted out of the active pool.</li>
    <li><strong>Maintenance-vs-resistance trade-off</strong> — Beulig 2025's headline
        finding: at high cell density, increased maintenance burden re-allocates
        proteome away from resistance functions toward maintenance, pushing cells
        into a persistence state. iML1515 can't see this; OxidizeME can.</li>
  </ul>
  <p>Operationally, the wrapper's extra outputs only populate under the real
  bridge:</p>
  <ul>
    <li><code>proteome_allocations: map[float]</code> — fraction per pathway
        (Σ MW × v_translation / Σ MW × v_translation_all)</li>
    <li><code>damaged_proteome_fraction: float</code> — Σ MW × v_ComplexFormation_damage
        / Σ MW × v_translation</li>
  </ul>
  <p>Both formulas come straight from Beulig 2025's Methods section
  (<a href="https://doi.org/10.1128/msystems.00323-25" target="_blank">doi:10.1128/msystems.00323-25</a>).
  See <code>OxidizeMEStep._compute_proteome_allocations</code> and
  <code>_compute_damaged_proteome</code> in <code>pbg_oxidizeme/processes.py</code>.</p>
</section>
"""


def render_install_obstacle() -> str:
    return """
<section class="card" style="border-left-color:#6b7280">
  <header><h3>Real bridge — install obstacle</h3></header>
  <p>The upstream OxidizeME stack runs on Python 2.7 / cobrapy 0.5.11 and requires
  <strong>qMINOS</strong> — a 128-bit-precision NLP solver — obtained by contacting
  <strong>Prof. M. A. Saunders, Stanford University</strong>. SoPlex (80-bit, academic-free)
  is an alternative. Both are required because ME-models are inherently ill-scaled.</p>
  <p>Install paths (in order of friction):</p>
  <ol>
    <li>Pull SBRG's Docker image:
        <code>docker pull sbrg/cobrame</code> — pre-installed full stack.</li>
    <li>Source build: clone SBRG/cobrame + SBRG/ecolime + SBRG/oxidizeme + SBRG/solvemepy
        into a Python 2.7 venv with cobrapy 0.5.11.</li>
    <li>Request qMINOS (or use SoPlex).</li>
  </ol>
  <p>Once installed, supply a pickled <code>StressME</code> model:</p>
  <pre>step = OxidizeMEStep(config={"me_model_pickle_path": "/path/to/stress_me.pickle", "strict": True})
sol = step.update({"external_concentrations": {"GLC[p]": 22.2, "OXYGEN-MOLECULE[p]": 0.21}})
print(sol["mu"], sol["proteome_allocations"], sol["damaged_proteome_fraction"])</pre>
</section>
"""


# ─── Top-level HTML template ────────────────────────────────────────────────


HTML_TEMPLATE = """<!DOCTYPE html><html lang="en"><head>
<meta charset="UTF-8"><title>pbg-oxidizeme — demo report</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
:root {{ --fg:#1f2937; --bg:#f9fafb; --card:#ffffff; --muted:#6b7280; --border:#e5e7eb; --accent:#6366f1; }}
* {{ box-sizing: border-box; }}
body {{ margin: 0; font: 14px/1.55 -apple-system, BlinkMacSystemFont, 'Segoe UI', system-ui, sans-serif;
       color: var(--fg); background: var(--bg); }}
header.page {{ position: sticky; top: 0; background: rgba(255,255,255,0.96); backdrop-filter: blur(10px);
              border-bottom: 1px solid var(--border); padding: 14px 28px; z-index: 50; }}
header.page h1 {{ margin: 0 0 4px; font-size: 20px; }}
header.page .meta {{ color: var(--muted); font-size: 12.5px; }}
nav.toc {{ display: flex; gap: 14px; font-size: 13px; margin-top: 8px; flex-wrap: wrap; }}
nav.toc a {{ color: #4338ca; text-decoration: none; }}
nav.toc a:hover {{ text-decoration: underline; }}
main {{ max-width: 1200px; margin: 28px auto; padding: 0 28px; }}
.banner {{ background:#dbeafe; border:1px solid #bfdbfe; color:#1e40af; border-radius:8px;
          padding:14px 18px; margin-bottom:22px; }}
.banner.warn {{ background:#fef3c7; border-color:#fde68a; color:#92400e; }}
.banner code {{ background: rgba(0,0,0,0.06); padding:1px 5px; border-radius:3px; }}
.card {{ background: var(--card); border: 1px solid var(--border); border-left: 5px solid var(--accent);
        border-radius: 10px; padding: 20px 24px; margin-bottom: 22px; }}
.card.warn {{ border-left-color: #ef4444; background:#fef2f2; }}
.card h2 {{ margin: 0 0 4px; font-size: 18px; }}
.card h3 {{ margin: 0 0 4px; font-size: 17px; }}
.card h4 {{ margin: 14px 0 8px; font-size: 14px; color:#374151; }}
.card .subtitle {{ color: var(--muted); margin: 0 0 8px; font-size: 13px; }}
.badges {{ display: flex; gap: 6px; flex-wrap: wrap; margin: 8px 0 12px; }}
.badge {{ display: inline-block; padding: 3px 10px; border-radius: 999px; font-size: 11.5px;
         background: #e5e7eb; color: #374151; font-weight: 600; }}
.badge.ok {{ background:#d1fae5; color:#065f46; }}
.badge.muted {{ background:#f3f4f6; color:#6b7280; font-weight:500; font-family: ui-monospace, monospace; }}
table {{ border-collapse: collapse; width: 100%; margin: 6px 0; font-size: 12.5px; }}
td, th {{ padding: 5px 10px; text-align: left; border-bottom: 1px solid var(--border); }}
.flux-up {{ color:#047857; font-family: ui-monospace, monospace; }}
.flux-down {{ color:#b91c1c; font-family: ui-monospace, monospace; }}
.two-col {{ display: grid; grid-template-columns: 1fr 1fr; gap: 18px; }}
.grid-2 {{ display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }}
.plotly {{ background: #fff; border-radius: 6px; }}
pre {{ background: #f3f4f6; padding: 10px 14px; border-radius: 6px; overflow: auto;
      font-size: 12px; line-height: 1.5; }}
code {{ background: rgba(0,0,0,0.06); padding: 1px 5px; border-radius: 3px;
        font-family: ui-monospace, Menlo, monospace; font-size: 12.5px; }}
.diagram {{ background: var(--card); border: 1px solid var(--border); border-radius: 10px;
           padding: 18px; text-align: center; margin-bottom: 22px; }}
.diagram img {{ max-width: 100%; height: auto; }}
ul, ol {{ padding-left: 22px; }}
ul li, ol li {{ margin: 4px 0; }}
@media (max-width: 800px) {{
  .two-col, .grid-2 {{ grid-template-columns: 1fr; }}
}}
</style>
</head><body>

<header class="page">
  <h1>pbg-oxidizeme — demo report</h1>
  <div class="meta">{generated_at} · process-bigraph wrapper for the OxidizeME ME-model (Yang et al. 2019 / Palsson lab)</div>
  <nav class="toc">
    <a href="#arch">Architecture</a>
    <a href="#iml1515">iML1515 simulation</a>
    <a href="#beulig">Beulig measured</a>
    <a href="#compare">iML1515 vs Beulig</a>
    <a href="#what-oxidizeme-adds">What OxidizeME adds</a>
    <a href="#install">Install obstacle</a>
  </nav>
</header>

<main>

<div class="banner">
  <strong>About this demo.</strong> Running <strong>now</strong>:
  <code>iML1515</code> (the M-model substrate that OxidizeME extends) via modern Python 3 cobrapy
  — three operating conditions, real FBA flux distributions —
  plus <strong>{n_reactors}</strong> reactors of real Beulig 2025 measured trajectories
  loaded from the digitized CSVs. <strong>Not running</strong>: OxidizeME itself (Python 2.7 / qMINOS
  install obstacle — see "Real bridge — install obstacle" below). When that stack
  is installed, the wrapper's port surface stays identical; the outputs gain
  <code>proteome_allocations</code> and <code>damaged_proteome_fraction</code>.
</div>

<section id="arch" class="diagram">
  <h2>Composite architecture</h2>
  {arch_img}
  <p style="text-align:left; max-width:760px; margin: 14px auto 0; color:#374151;">
    The <code>OxidizeMEStep</code> consumes <code>external_concentrations</code> (EcoCyc-keyed,
    mM), <code>ros_concentrations</code> + <code>metal_concentrations</code> (nmol/L), and an optional
    <code>mu_fixed</code>. It emits <code>mu</code>, <code>external_exchange_fluxes</code>,
    <code>proteome_allocations</code>, <code>damaged_proteome_fraction</code>, and
    <code>solver_status</code>. Ports conform to the cell-side interface contract used by
    v2ecoli's multiscale-bioprocess investigation — drop-in interchangeable with
    v2ecoli baseline under the same reactor composite.
  </p>
</section>

<section id="iml1515">
<h2 style="margin:24px 0 12px;">iML1515 simulation (live)</h2>
<p style="color:#374151;">
  iML1515 is the genome-scale M-model OxidizeME is built on. Running it on its own (no expression
  layer) at three operating conditions previews where the wrapper's flux outputs come from — the
  expression / ROS layer modulates these but does not redo the metabolic flux balance.
</p>
{iml1515_cards}
</section>

<section id="beulig">
{beulig_section}
</section>

<section id="compare">
{compare_section}
</section>

<section id="what-oxidizeme-adds">
{what_adds}
</section>

<section id="install">
{install_obstacle}
</section>

</main>
</body></html>
"""


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("[demo] running iML1515 configs…")
    iml_results = run_iml1515_configs()
    print(f"[demo] iML1515: {len(iml_results)} configs, μ = {[r['mu'] for r in iml_results]}")

    print("[demo] loading Beulig trajectories…")
    beulig = load_beulig_trajectories()
    if beulig["available"]:
        print(f"[demo] Beulig: {len(beulig['reactor_ids'])} reactors")
    else:
        print(f"[demo] Beulig: data not available at {beulig.get('missing_path')}")

    print("[demo] rendering bigraph PNG…")
    img_uri = render_bigraph_png()
    arch_img = (
        f'<img src="{img_uri}" alt="composite architecture">'
        if img_uri else "<p>(bigraph-viz unavailable)</p>"
    )

    iml1515_cards = "\n".join(render_iml1515_card(r) for r in iml_results)
    beulig_section = render_beulig_section(beulig)
    compare_section = render_iml1515_vs_beulig(iml_results, beulig)
    what_adds = render_what_oxidizeme_adds()
    install_obstacle = render_install_obstacle()

    html = HTML_TEMPLATE.format(
        generated_at=time.strftime("%Y-%m-%d %H:%M:%S %Z"),
        n_reactors=len(beulig.get("reactor_ids", [])),
        arch_img=arch_img,
        iml1515_cards=iml1515_cards,
        beulig_section=beulig_section,
        compare_section=compare_section,
        what_adds=what_adds,
        install_obstacle=install_obstacle,
    )
    REPORT_PATH.write_text(html, encoding="utf-8")
    print(f"[demo] wrote {REPORT_PATH} ({len(html)/1024:.0f} KB)")
    webbrowser.open("file://" + str(REPORT_PATH.resolve()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
