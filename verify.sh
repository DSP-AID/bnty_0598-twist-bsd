#!/usr/bin/env bash
# Reproduce the verification from a clean clone (needs python3 + mpmath, PARI/GP optional).
set -e
cd "$(dirname "$0")"
python3 -c "import mpmath" 2>/dev/null || pip install --quiet mpmath

echo "==> [1/2] engine self-test (8 analytic anchors incl. L(E,1)=Gamma(1/4)^2/(4*sqrt(pi))):"
python3 src/engine.py

echo
echo "==> [2/2] full cross-validation of the 404-curve dataset (root numbers, ranks,"
echo "    parity, ball-vs-PARI agreement, BSD quotient = integer square):"
python3 src/merge.py

echo
echo "ALL CHECKS PASSED"
