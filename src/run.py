"""Certified analytic BSD-side data for E^(d): y^2 = x^3 - d^2 x.

Usage: python run.py D OUT.jsonl   (odd squarefree 1 <= d <= D)

For each d emits one JSON line: conductor N=32d^2, root number w=(2/d)-type
formula sign, certified ball [val-rad, val+rad] for L(E^(d),1) (w=+1) or
L'(E^(d),1) (w=-1), the truncation index M, and the analytic rank verdict
r_an = 0 (w=+1, ball excludes 0) or 1 (w=-1, ball excludes 0).

Run on Linux server; cross-checked against PARI/GP lfun/ellrootno (suite.gp).
"""

import json
import os
import sys

import engine
import mpmath
from engine import (
    L_value,
    _M_cut,
    kronecker_m2,
    squarefree_odd_upto,
    twisted_coeffs,
)

# MSCALE: multiply truncation index (and working precision) for hard curves
# whose |L(1)| falls below the default ball radius (~1e-38).
MSCALE = int(os.environ.get("MSCALE", "1"))
if MSCALE > 1:
    mpmath.mp.dps = 45 * MSCALE
    engine._CUT_SCALE = MSCALE


def cert(d):
    N = 32 * d * d
    w = kronecker_m2(d)
    M = _M_cut(N) * MSCALE + 5
    an = twisted_coeffs(d, M)
    if w == 1:
        v, r, M = L_value(N, an, w)
    else:
        v, r, M = L_value(N, an, w, deriv=True)
    nonzero = abs(v) > r
    return {
        "d": d,
        "N": N,
        "w": w,
        "val": mpmath.nstr(v, 40),
        "rad": mpmath.nstr(r, 3),
        "M": M,
        "ball_excl_0": bool(nonzero),
        "r_an": (0 if w == 1 else 1) if nonzero else None,
    }


if __name__ == "__main__":
    D = int(sys.argv[1])
    ds = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else squarefree_odd_upto(D)
    with open(sys.argv[2], "w", encoding="utf-8") as out:
        for d in ds:
            rec = cert(d)
            out.write(json.dumps(rec) + "\n")
            out.flush()
            print(f"d={d} w={rec['w']:+d} r_an={rec['r_an']}", file=sys.stderr)
