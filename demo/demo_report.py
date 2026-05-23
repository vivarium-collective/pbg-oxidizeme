"""Demo report for pbg-oxidizeme.

Runs the steady-state composite in non-strict mode (so it works without the
upstream OxidizeME solver installed) and renders a self-contained HTML
report. When the upstream stack IS installed and OXIDIZEME_PICKLE points at
a pickled StressME model, the same script exercises the REAL bridge end-to-
end and the report carries the real solution.

Output: demo/report.html
"""

from __future__ import annotations

import base64
import json
import os
import sys
import time
import webbrowser
from pathlib import Path
from typing import Any

# Allow running directly from the repo root
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = REPO_ROOT / "demo"
REPORT_PATH = OUT_DIR / "report.html"


CONFIGS: list[dict[str, Any]] = [
    {
        "id": "aerobic-baseline",
        "title": "Aerobic batch, baseline ROS",
        "subtitle": "WT BW25113 at M9-glucose, near-saturation O2, textbook ROS",
        "description": (
            "The reference run. OxidizeME is solved at the textbook M9-glucose "
            "exchange bounds with paper-default H2O2 / O2- intracellular "
            "concentrations. When the real solver is installed, this gives the "
            "WT μ_max baseline against which other configurations are compared."
        ),
        "params": {
            "glucose_mM": 22.2, "oxygen_mM": 0.21, "ammonium_mM": 19.0,
            "h2o2_c_nmolL": 0.2, "o2s_c_nmolL": 0.05,
        },
        "accent": "#6366f1",
    },
    {
        "id": "low-o2",
        "title": "O2-limited (mimics Beulig fed-batch onset)",
        "subtitle": "Low dissolved O2 (~5% of saturation)",
        "description": (
            "Reactor-limited dissolved O2 — the regime where Beulig 2025 sees "
            "overflow + transcriptional stress transitions. OxidizeME should "
            "predict reduced μ + shift toward fermentation byproducts."
        ),
        "params": {
            "glucose_mM": 22.2, "oxygen_mM": 0.01, "ammonium_mM": 19.0,
            "h2o2_c_nmolL": 0.2, "o2s_c_nmolL": 0.05,
        },
        "accent": "#8b5cf6",
    },
    {
        "id": "elevated-ros",
        "title": "Elevated ROS (paraquat-like challenge)",
        "subtitle": "10× baseline H2O2 + O2- to invoke damage/repair",
        "description": (
            "Pushes the ROS-damage subnet OxidizeME is designed for. Damaged "
            "proteome fraction should rise; proteome allocation should shift "
            "toward damage-repair complexes; μ_max may drop."
        ),
        "params": {
            "glucose_mM": 22.2, "oxygen_mM": 0.21, "ammonium_mM": 19.0,
            "h2o2_c_nmolL": 2.0, "o2s_c_nmolL": 0.5,
        },
        "accent": "#ec4899",
    },
]


def run_config(cfg: dict[str, Any]) -> dict[str, Any]:
    """Build and run a single composite. Returns a snapshot record."""
    from pbg_oxidizeme.composites import oxidizeme_steady_state
    from process_bigraph import Composite, allocate_core

    core = allocate_core()
    # Honour OXIDIZEME_PICKLE env var → real-bridge solve when present.
    pickle_path = os.environ.get("OXIDIZEME_PICKLE", "")
    strict = bool(pickle_path)   # if pickle is set, run the real solver

    doc = oxidizeme_steady_state(
        core=core,
        me_model_pickle_path=pickle_path,
        strict=strict,
        **cfg["params"],
    )
    sim = Composite({"state": doc}, core=core)

    t0 = time.perf_counter()
    sim.run(1.0)   # one tick — the Step computes the solve
    wall = time.perf_counter() - t0

    # Read state out of the Composite via gather_emitter_results
    from process_bigraph.emitter import gather_emitter_results
    emitted = gather_emitter_results(sim) or {}
    # Pick the most recent snapshot from any emitter path
    last = {}
    for _path, snapshots in emitted.items():
        if snapshots:
            last = snapshots[-1]
    stores = last
    return {
        "id": cfg["id"],
        "title": cfg["title"],
        "subtitle": cfg["subtitle"],
        "description": cfg["description"],
        "params": cfg["params"],
        "wall_s": wall,
        "accent": cfg["accent"],
        "solver_status": stores.get("solver_status", "unknown"),
        "mu": stores.get("mu", 0.0),
        "exchange_fluxes": stores.get("exchange_fluxes", {}),
        "proteome_allocations": stores.get("proteome_allocations", {}),
        "damaged_proteome": stores.get("damaged_proteome", 0.0),
    }


def render_bigraph_png() -> str | None:
    """Render a bigraph-viz PNG of the steady-state composite. Returns a
    data: URI, or None if bigraph-viz isn't installed.
    """
    try:
        from bigraph_viz import plot_bigraph   # type: ignore
    except ImportError:
        return None
    from pbg_oxidizeme.composites import oxidizeme_steady_state

    doc = oxidizeme_steady_state(core=None, strict=False)
    out_dir = OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        plot_bigraph(
            state=doc,
            out_dir=str(out_dir),
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
    png = out_dir / "bigraph.png"
    if not png.exists():
        return None
    with open(png, "rb") as fh:
        return "data:image/png;base64," + base64.b64encode(fh.read()).decode()


def html_card(record: dict[str, Any]) -> str:
    mu = record["mu"]
    mu_s = "NaN" if (isinstance(mu, float) and mu != mu) else f"{mu:.3f}"
    status = record["solver_status"]
    if status == "upstream_missing":
        status_badge = '<span class="badge badge-empty">upstream not installed</span>'
    elif status == "optimal":
        status_badge = '<span class="badge badge-ok">optimal</span>'
    else:
        status_badge = f'<span class="badge badge-warn">{status}</span>'

    flux_rows = "".join(
        f"<tr><td>{k}</td><td>{v:.4g}</td></tr>"
        for k, v in record["exchange_fluxes"].items()
    ) or '<tr><td colspan="2" class="muted">no fluxes (upstream not installed)</td></tr>'

    pa_rows = "".join(
        f"<tr><td>{k}</td><td>{v:.3%}</td></tr>"
        for k, v in record["proteome_allocations"].items()
    ) or '<tr><td colspan="2" class="muted">no allocations (upstream not installed)</td></tr>'

    params_json = json.dumps(record["params"], indent=2)

    return f"""
<section class="card" style="border-left-color:{record['accent']}">
  <header>
    <h2>{record['title']}</h2>
    <p class="subtitle">{record['subtitle']}</p>
    <div class="badges">
      {status_badge}
      <span class="badge">μ = {mu_s} 1/h</span>
      <span class="badge">wall = {record['wall_s']*1000:.1f} ms</span>
    </div>
  </header>
  <p>{record['description']}</p>
  <details><summary>Input parameters</summary><pre>{params_json}</pre></details>
  <h3>External exchange fluxes</h3>
  <table><thead><tr><th>EcoCyc ID</th><th>mmol/(gDW·h)</th></tr></thead><tbody>{flux_rows}</tbody></table>
  <h3>Proteome allocations</h3>
  <table><thead><tr><th>Pathway</th><th>Fraction</th></tr></thead><tbody>{pa_rows}</tbody></table>
</section>
"""


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"><title>pbg-oxidizeme demo</title>
<style>
:root {{ --fg:#1f2937; --bg:#f9fafb; --card:#ffffff; --muted:#6b7280; --border:#e5e7eb; }}
* {{ box-sizing: border-box; }}
body {{ margin: 0; font: 14px/1.5 -apple-system, BlinkMacSystemFont, 'Segoe UI', system-ui, sans-serif;
       color: var(--fg); background: var(--bg); }}
header.page {{ position: sticky; top: 0; background: rgba(255,255,255,0.95); backdrop-filter: blur(8px);
              border-bottom: 1px solid var(--border); padding: 12px 24px; z-index: 10; }}
header.page h1 {{ margin: 0; font-size: 18px; }}
header.page .subtitle {{ color: var(--muted); font-size: 13px; }}
main {{ max-width: 1100px; margin: 24px auto; padding: 0 24px; }}
.banner {{ background: #fef3c7; border: 1px solid #fde68a; border-radius: 6px; padding: 12px 16px;
          margin-bottom: 18px; }}
.banner.ok {{ background:#d1fae5; border-color:#a7f3d0; }}
.card {{ background: var(--card); border: 1px solid var(--border); border-left: 4px solid var(--muted);
        border-radius: 8px; padding: 18px 22px; margin-bottom: 18px; }}
.card h2 {{ margin: 0 0 4px; font-size: 18px; }}
.card .subtitle {{ color: var(--muted); margin: 0 0 8px; font-size: 13px; }}
.badges {{ display: flex; gap: 8px; flex-wrap: wrap; margin: 6px 0 8px; }}
.badge {{ display: inline-block; padding: 2px 8px; border-radius: 999px; font-size: 11px;
         background: #e5e7eb; color: #374151; font-weight: 600; }}
.badge-ok {{ background:#d1fae5; color:#065f46; }}
.badge-warn {{ background:#fee2e2; color:#991b1b; }}
.badge-empty {{ background:#e0e7ff; color:#3730a3; }}
.muted {{ color: var(--muted); font-style: italic; }}
table {{ border-collapse: collapse; width: 100%; margin: 8px 0 16px; font-size: 13px; }}
td, th {{ padding: 4px 8px; text-align: left; border-bottom: 1px solid var(--border); }}
th {{ background: #f3f4f6; }}
details {{ margin: 8px 0; }}
summary {{ cursor: pointer; color: #4f46e5; font-size: 13px; }}
pre {{ background: #f3f4f6; padding: 8px 12px; border-radius: 6px; overflow: auto; font-size: 12px; }}
.diagram {{ background: var(--card); border: 1px solid var(--border); border-radius: 8px; padding: 18px;
           text-align: center; margin-bottom: 18px; }}
.diagram img {{ max-width: 100%; height: auto; }}
.diagram h2 {{ margin: 0 0 12px; text-align: left; font-size: 16px; }}
</style>
</head><body>

<header class="page">
  <h1>pbg-oxidizeme — demo report</h1>
  <div class="subtitle">{generated_at} · process-bigraph wrapper for the OxidizeME ME-model</div>
</header>

<main>
{banner}
{diagram}
{cards}

<section class="card">
<h2>Cell-side interface contract</h2>
<p>This Step satisfies the substitutability contract described at
<code>references/expert/cell_side_interface_contract.md</code> in the v2ecoli
multiscale-bioprocess investigation. Inputs and outputs use EcoCyc-style
compartment-tagged IDs (e.g. <code>GLC[p]</code>) and bare composable types
where applicable, so this engine can be swapped in for v2ecoli baseline
under the same reactor composite without re-wiring.</p>
</section>

<section class="card">
<h2>Real bridge — install obstacle</h2>
<p>The upstream OxidizeME stack is Python 2.7 / cobrapy 0.5.11 era; qMINOS
is proprietary (request from Prof. M. A. Saunders, Stanford). The wrapper
runs the real solver when the stack is installed; otherwise it returns a
labelled empty solution so dashboard wiring still works. See the README
for the install path (SBRG Docker image / qMINOS request / SoPlex
alternative).</p>
</section>
</main>
</body></html>
"""


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    results = [run_config(cfg) for cfg in CONFIGS]

    upstream_present = any(r["solver_status"] not in {"upstream_missing", "error:OxidizeMEMissingError"} for r in results)
    if upstream_present:
        banner = '<div class="banner ok">✓ Upstream OxidizeME stack detected — real-bridge solutions below.</div>'
    else:
        banner = (
            '<div class="banner">⚠ Upstream OxidizeME stack not installed — '
            'solutions below are labelled empty (NaN scalars, empty dicts). '
            'This page exercises the wrapper\'s port surface but NOT the science. '
            'See README "Real bridge — install obstacle".</div>'
        )

    img_uri = render_bigraph_png()
    if img_uri:
        diagram = f'<div class="diagram"><h2>Composite architecture</h2><img src="{img_uri}" alt="bigraph"></div>'
    else:
        diagram = ""

    cards = "\n".join(html_card(r) for r in results)

    html = HTML_TEMPLATE.format(
        generated_at=time.strftime("%Y-%m-%d %H:%M:%S %Z"),
        banner=banner,
        diagram=diagram,
        cards=cards,
    )
    REPORT_PATH.write_text(html, encoding="utf-8")
    print(f"wrote {REPORT_PATH}")
    webbrowser.open("file://" + str(REPORT_PATH.resolve()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
