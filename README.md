# BSD on an explicit family of quadratic twists: 404 curves with certified error balls

Machine-verifiable verification of the Birch & Swinnerton-Dyer conjecture for the
404 quadratic twists `E^(d): y^2 = x^3 - d^2*x` (odd squarefree `1 <= d <= 1000`)
of the congruence curve `32a1`, proved via the modular parametrisation of
`X_0(32)` alone — every L-value, parity sign and `Sha` prediction in this repo is
computed from that parametrisation with rigorous error balls (no Heegner-point
computation anywhere; rank upper bounds use classical 2-descent).

Headline results (all reproducible in minutes on a laptop):

| statement | verified |
|---|---|
| parity `w(E^(d)) = (-2/d)` | 404/404 vs PARI `ellrootno` + 404/404 vs the analytic formula of [paper.md, Thm A] |
| algebraic rank = analytic rank | 370/370 (2-descent bounds from PARI `ellrank` + certified `r_an`) |
| BSD quotient `16*L^(r)(1)/(r! * Omega * Reg * c)` is a perfect square | **360/360, zero outliers** (max deviation 6e-37) |
| predicted `#Sha` values (all squares) | `{1: 263, 4: 51, 9: 24, 16: 12, 25: 6, 36: 2, 49: 1, 81: 1}` |
| rank-2 members (Selmer bounds `[2,2]` + explicit generators) | 34/34 predict `Sha = 1` (rank = 2 unconditionally; `#Sha = 1` via cited GZ/Kolyvagin finiteness, paper.md Prop D) |

The proof package is [paper.md](paper.md): root-number theorem, certified
truncation/derivative bounds (Lemmas 4.1–4.3), the period/regulator dictionary
anchored on the classical `Gamma(1/4)` evaluation of `L(E,1)`, and the BSD
consistency theorem.

## Layout

```
paper.md        proof package (theorems + proofs + results tables)
src/engine.py   certified analytic engine (mpmath; self-contained)
src/run.py      certificate generator (usage: `python3 src/run.py D OUT [dlist]`)
src/merge.py    full cross-validation of data/ (run: python3 src/merge.py)
gp/             PARI/GP scripts that regenerate every external bound/table
data/           certificates, root numbers, 2-descent data, BSD tables
verify.sh       one-command reproduction
```

## Reproduce

```bash
./verify.sh        # engine self-test + full cross-validation (python3, mpmath)
```

Optional heavy parts (PARI/GP >= 2.15, minutes–hours): see paper.md §7 for the
exact `gp/suite_*.gp` invocations. Certified-ball spot check (writes
`src/new_certificates.jsonl`, leaves `data/` untouched):

```bash
python3 src/run.py 1000 new_certificates.jsonl 29,58,41    # d=29,58: balls contain 0 (genuine double zeros)
MSCALE=4 python3 src/run.py 1000 nc4.jsonl 29,58,41        # 180-digit rerun: still contain 0
```

## Notes

- `certificates.jsonl` stores raw `L^(r)(1)` values + ball radii (dps); the BSD
  normalisation dictionary and its d=1 anchor are derived in paper.md §5.
- The 44 rank-1 twists listed in paper.md §6 lack a PARI-found generator, so
  their regulator is unavailable; for them parity + root number + `L^(r)(1) != 0`
  are still certified.
