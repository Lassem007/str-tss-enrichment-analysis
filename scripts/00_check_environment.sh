#!/usr/bin/env bash
set -euo pipefail

printf "Checking required command-line tools...\n"
for tool in bedtools awk sort cut grep wc gunzip wget python3; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "ERROR: $tool is not installed or not in PATH" >&2
    exit 1
  fi
  echo "OK: $tool"
done

printf "\nBEDTools version:\n"
bedtools --version

printf "\nPython packages required for visualization:\n"
python3 - <<'PY'
import importlib
for pkg in ["numpy", "pandas", "matplotlib"]:
    try:
        importlib.import_module(pkg)
        print(f"OK: {pkg}")
    except ImportError:
        print(f"MISSING: {pkg}. Install with: pip3 install pandas numpy matplotlib")
PY
