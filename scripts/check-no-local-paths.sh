#!/usr/bin/env bash
# Guard against local laptop paths leaking into committed config.
#
# A recurring failure mode in pbg workspaces: a local editable install such as
#   vivarium-dashboard = { path = "/Users/<you>/code/vivarium-dashboard", ... }
# gets committed to pyproject.toml. uv then synthesizes a `file:///Users/...`
# URL that exists only on the author's laptop, so CI (and any fresh clone)
# fails with "Distribution not found at: file:///Users/...".
#
# This script greps every TRACKED file for absolute local paths and fails fast
# with a clear message. It excludes itself (it necessarily contains the
# patterns it searches for).
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

# Patterns that should never appear in committed files. Assembled at runtime so
# this script does not match itself when it scans the tree.
patterns=(
  "file:""///"        # local file:// URLs (uv expands path= sources to these)
  "/Users""/"         # macOS home dirs
  "/home""/"          # linux home dirs (laptop checkouts)
)

# Files allowed to mention these patterns (docs that show them as examples).
exclude_pathspecs=(
  ":!scripts/check-no-local-paths.sh"
  ":!.github/workflows/workspace-ci.yml"
)

found=0
for pat in "${patterns[@]}"; do
  if git grep -n -- "$pat" -- "${exclude_pathspecs[@]}" 2>/dev/null; then
    found=1
  fi
done

if [ "$found" -ne 0 ]; then
  echo ""
  echo "ERROR: a local absolute path leaked into a tracked file (see above)."
  echo "       Replace local 'path = \"/Users/...\"' deps with a git source, e.g.:"
  echo "         vivarium-dashboard = { git = \"https://github.com/vivarium-collective/vivarium-dashboard.git\", branch = \"main\" }"
  echo "       For local development use an editable install in your venv instead:"
  echo "         uv pip install -e ../vivarium-dashboard"
  exit 1
fi

echo "OK: no local absolute paths in tracked files."
