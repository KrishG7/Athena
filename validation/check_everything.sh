#!/usr/bin/env bash
# Run every check this project has, in dependency order, and say plainly which
# passed. Offline checks always run. Live checks need FalkorDB and network, and
# are skipped (not failed) when those are missing.
#
#   bash validation/check_everything.sh            # offline checks only
#   bash validation/check_everything.sh --live     # also run live/integration
#
# Override the interpreter with PYTHON=/path/to/python.

set -uo pipefail
cd "$(dirname "$0")/.."
REPO="$PWD"
LIVE=0
[ "${1:-}" = "--live" ] && LIVE=1

# --- locate an interpreter that has the dependencies -------------------------
if [ -n "${PYTHON:-}" ]; then
    :
elif [ -x "$REPO/.venv/bin/python" ]; then
    PYTHON="$REPO/.venv/bin/python"
elif [ -x "$REPO/../Athena_tier_2/.venv/bin/python" ]; then
    PYTHON="$REPO/../Athena_tier_2/.venv/bin/python"
else
    PYTHON="$(command -v python3 || true)"
fi
[ -x "$PYTHON" ] || { echo "no usable python found; set PYTHON=/path/to/python"; exit 2; }

PASS=0; FAIL=0; SKIP=0
ok()   { printf '  \033[32mPASS\033[0m  %s\n' "$1"; PASS=$((PASS+1)); }
bad()  { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; FAIL=$((FAIL+1)); }
skip() { printf '  \033[33mSKIP\033[0m  %s\n' "$1"; SKIP=$((SKIP+1)); }
hdr()  { printf '\n\033[1m%s\033[0m\n' "$1"; }

echo "Athena — full check"
echo "repo:   $REPO"
echo "python: $PYTHON"

# --- 1. the code that will actually be imported ------------------------------
hdr "1. Import resolution"
# secondlook is a PEP 420 namespace package. If an older checkout is on the
# path (an editable install pointing elsewhere), modules can silently resolve
# to it. Confirm they come from THIS repo before trusting anything below.
RESOLVED=$("$PYTHON" - <<'PY' 2>/dev/null
import sys; sys.path.insert(0, "src")
import secondlook.pipeline as p
print(p.__file__)
PY
)
case "$RESOLVED" in
    "$REPO"/src/secondlook/pipeline.py) ok "secondlook resolves to this repo" ;;
    "") bad "could not import secondlook.pipeline — dependencies missing?" ;;
    *)  bad "secondlook resolves to ANOTHER checkout: $RESOLVED" ;;
esac

# --- 2. unit tests -----------------------------------------------------------
hdr "2. Unit tests (offline)"
if OUT=$("$PYTHON" -m pytest -q 2>&1); then
    ok "$(echo "$OUT" | tail -1 | sed 's/\x1b\[[0-9;]*m//g')"
else
    bad "unit tests failed"; echo "$OUT" | tail -20
fi

# --- 3. self-containment -----------------------------------------------------
hdr "3. Self-containment (would a fresh clone work?)"
# Re-run the suite with any sibling checkout removed from sys.path. A module
# missing from this repo would otherwise import from the sibling and the suite
# would pass on code a fresh clone does not contain.
if OUT=$("$PYTHON" - -q <<'PY' 2>&1
import sys
from pathlib import Path
here = str(Path.cwd() / "src")
sys.path[:] = [p for p in sys.path
               if not (p.rstrip("/").endswith("/src") and p != here)]
import pytest
sys.exit(pytest.main(sys.argv[1:]))
PY
); then
    ok "suite passes with sibling checkouts off sys.path — $(echo "$OUT" | tail -1 | sed 's/\x1b\[[0-9;]*m//g')"
else
    bad "suite depends on a source tree outside this repo"; echo "$OUT" | tail -20
fi

# --- 4. proximity verification ----------------------------------------------
hdr "4. Proximity measurement — pre-registered criteria C1-C4"
# Independent re-derivation from raw PDB coordinates. Does not import
# proximity.py: sharing code with the thing under test would prove nothing.
if [ -d validation/cache ] && ls validation/cache/*.json >/dev/null 2>&1; then
    if OUT=$("$PYTHON" validation/verify_proximity.py --cache validation/cache 2>&1); then
        echo "$OUT" | grep -E '^C[1-4] ' | sed 's/^/        /'
        ok "all four criteria pass"
    else
        echo "$OUT" | grep -E '^C[1-4] ' | sed 's/^/        /'
        bad "a pre-registered criterion failed — see docs/proximity-validation-plan.md"
    fi
else
    skip "no validation/cache — run: $PYTHON validation/run_gold_standard.py"
fi

# --- 5. gate triage ----------------------------------------------------------
hdr "5. Gate triage (which cases can produce a delta at all)"
# Seconds, not hours: the covalent gate, structure/protomer choice and the 8 A
# contact range all decide before docking starts, and they explain most of why a
# case does or does not produce a number. Reported, not pass/failed -- a refusal
# here is usually correct behaviour.
if OUT=$("$PYTHON" validation/triage.py 2>&1); then
    echo "$OUT" | sed 's/^/        /'
    ok "triage ran"
else
    echo "$OUT" | tail -10 | sed 's/^/        /'
    bad "triage failed"
fi

# --- 6. live checks ----------------------------------------------------------
hdr "6. Live checks (network + FalkorDB)"
if [ "$LIVE" -eq 0 ]; then
    skip "not requested — re-run with --live"
else
    if "$PYTHON" -c "
import socket,sys
s=socket.socket(); s.settimeout(2)
sys.exit(0 if s.connect_ex(('localhost',6379))==0 else 1)" 2>/dev/null; then
        ok "FalkorDB reachable on localhost:6379"
    else
        skip "FalkorDB not on localhost:6379 — 'docker compose up -d'; graph tests will skip"
    fi
    if OUT=$("$PYTHON" -m pytest -m integration -q 2>&1); then
        ok "integration tests — $(echo "$OUT" | tail -1 | sed 's/\x1b\[[0-9;]*m//g')"
    else
        bad "integration tests failed"; echo "$OUT" | tail -25
    fi
fi

hdr "Summary"
printf '  %d passed, %d failed, %d skipped\n' "$PASS" "$FAIL" "$SKIP"
if [ "$FAIL" -gt 0 ]; then
    echo "  Something is wrong. ISSUES.md maps every known problem to a fix."
    exit 1
fi
echo "  All checks that ran passed."
echo
echo "  Not covered here (needs a long docking run, see ISSUES.md §15):"
echo "    - docking noise floor at exhaustiveness 32"
echo "    - gold-standard re-run since the reproducibility fix"
