# Contributing to pbg-oxidizeme

## Development setup

uv is required. Install with `brew install uv` or `pip install uv`.

```bash
uv venv .venv
source .venv/bin/activate
uv pip install -e ".[dev]"
pytest
```

The `bigraph-schema` + `process-bigraph` deps install here. The OxidizeME
stack (cobrame, ecolime, qminospy, oxidizeme) is environmental — see the
README for the install path. Tests that exercise the real solver skip
gracefully when the upstream stack is absent.

## Adding a config knob

OxidizeMEStep's `config_schema` is the build-time surface; per-call state
lives on input ports. When adding a config knob, ask:

- Could a sibling Process sensibly write this each step?  → input port.
- Is it a calibration / model-loading setting fixed at construction? → config.

See [Port Design](https://github.com/vivarium-collective/pbg-superpowers/blob/main/skills/pbg-expert/SKILL.md#port-design)
in the pbg-expert skill for the rationale.

## Releasing to PyPI

Tag a commit with `git tag v<VERSION>` and push the tag. The
`.github/workflows/release.yml` workflow publishes to PyPI automatically
using trusted publishing (no tokens needed after initial setup).

PyPI trusted publishing must be configured once per repo — see
https://docs.pypi.org/trusted-publishers/ and the pbg-superpowers
distribution guide.
