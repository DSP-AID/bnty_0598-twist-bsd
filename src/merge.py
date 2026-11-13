"""merge.py - cross-validate engine certificates vs PARI/GP ground truth (bnty_0598).

Inputs (same dir):
  certificates.jsonl / certs_extra.jsonl  engine balls for L^{(r)}(E^(d),1)
  rootno.txt    "R d w"            ellrootno cross-check
  bsd_d1000.txt "B d N omega1 tam [2, 2] w V"  PARI conductor/period/tamagawa/sign/value
  rank_d250.txt / rank_d1000.txt  "K d rk [heights]"  ellrank + canonical heights
  reg.txt       "G d [rlo, rhi] Reg2 L2"  rank-2 regulators + L''(1) (optional)

Output: results.csv + summary stats. Pure stdlib (Decimal arithmetic).

BSD quotient convention (verified on d=1,3,5,7):
  #Sha_pred = 16 * L^{(r)}(1)/r! / (omega1 * Reg * tam),  tam = elltamagawa (incl c_inf=2)
  equivalent to BSD with Omega_BSD(E^(d)) = sqrt(d)*Omega_BSD(32a1).
"""
import csv
import json
import math
import os
from decimal import ROUND_FLOOR, Decimal, getcontext

getcontext().prec = 60
HERE = os.path.dirname(os.path.abspath(__file__))


def _path(name):
    for d in (HERE, os.path.join(HERE, os.pardir, "data")):
        p = os.path.join(d, name)
        if os.path.exists(p):
            return p
    return None


def read_lines(name):
    p = _path(name)
    if p is None:
        return []
    with open(p, encoding="utf-8") as f:
        return [ln.strip() for ln in f if ln.strip()]


def dec(s):
    return Decimal("".join(s.split()))  # strip PARI print's inner spaces (E-notation)


def kronecker_m2(d):
    return -1 if d % 8 in (5, 7) else 1


def load_certs():
    certs = {}
    for name in ["certificates.jsonl", "certs_extra.jsonl"]:
        for ln in read_lines(name):
            j = json.loads(ln)
            d = j["d"]
            if name == "certs_extra.jsonl" or d not in certs:  # extra supersedes; base first-wins
                certs[d] = j
    return certs


def main():
    certs = load_certs()

    rootno = {}
    for ln in read_lines("rootno.txt"):
        p = ln.split()
        if p[0] == "R" and int(p[1]) not in rootno:
            rootno[int(p[1])] = int(p[2])

    bsd = {}
    for ln in read_lines("bsd_d1000.txt"):
        p = ln.split()
        if p[0] != "B":
            continue
        d = int(p[1])
        if d not in bsd:
            bsd[d] = {"N": int(p[2]), "om1": dec(p[3]), "tam": Decimal(p[4]),
                      "w": int(p[7]), "V": dec(" ".join(p[8:]))}

    rank = {}
    for name in ["rank_d250.txt", "rank_d1000.txt"]:
        for ln in read_lines(name):
            p = ln.split()
            if p[0] != "K":
                continue
            d = int(p[1])
            rk = int(p[2])
            h1 = dec(p[3].strip("[], ")) if rk == 1 and len(p) > 3 and p[3].strip("[], ") else None
            rank[d] = {"rk": rk, "h1": h1}

    reg2 = {}
    for ln in read_lines("reg.txt"):
        p = ln.split()
        if p[0] != "G":
            continue
        # "G d;r;rbounds;Reg2;L2" - semicolon-delimited, immune to PARI space quirks
        f = "".join(ln.split()[1:]).split(";")
        reg2[int(f[0])] = {"r": int(f[1]), "rmin": int(f[2]),
                           "reg": dec(f[3]), "L2": dec(f[4])}

    ds = sorted(set(certs) & set(bsd))
    rows, st = [], {"n": 0, "Nbad": [], "wbad": [], "vdisagree": [], "ran": {0: 0, 1: 0, None: 0},
                    "rankagree": 0, "rankdis": [], "sqc": 0, "certsq": 0, "nonsq": [],
                    "sha": {}, "maxvd": Decimal(0), "r2list": [], "norank": [], "reg2miss": []}
    for d in ds:
        c, b = certs[d], bsd[d]
        row = {"d": d, "N": b["N"], "w": b["w"], "w_formula": kronecker_m2(d),
               "w_rootno": rootno.get(d, "")}
        if b["N"] != 32 * d * d:
            st["Nbad"].append(d)
        if not (b["w"] == kronecker_m2(d) == c["w"] == rootno.get(d, b["w"])):
            st["wbad"].append(d)
        val = dec(c["val"])
        rad = dec(c["rad"])
        vd = abs(val - b["V"])
        st["maxvd"] = max(st["maxvd"], min(vd, Decimal("1e-30")))
        if vd > rad + Decimal("1e-28") * max(abs(b["V"]), Decimal(1)):
            st["vdisagree"].append((d, str(vd)[:12]))
        decisive = abs(val) > rad
        r_an = c.get("r_an")
        if decisive and r_an is None:
            r_an = 0 if b["w"] == 1 else 1
        st["ran"][r_an] += 1

        rk = rank.get(d, {}).get("rk")
        row["pari_rank"] = rk if rk is not None else "?"
        if rk is None:
            st["norank"].append(d)
        elif r_an is not None and rk != r_an:
            st["rankdis"].append((d, r_an, rk))
        elif r_an is not None:
            st["rankagree"] += 1

        # BSD quotient per rank
        reg, Vbsd, fact = None, b["V"], 1
        if rk == 0:
            reg, mode = Decimal(1), "r0"
        elif rk == 1 and rank[d]["h1"]:
            reg, mode = rank[d]["h1"], "r1"
        elif rk == 2:
            st["r2list"].append(d)
            g = reg2.get(d)
            if not g or not g["L2"]:
                st["reg2miss"].append(d)
                rows.append(row)
                continue
            reg, mode = g["reg"], "r2"
            Vbsd, fact = g["L2"], 2  # L''(1)/2!
        else:
            if rk == 1:
                st.setdefault("nogen", []).append(d)
            if rk is not None and r_an is None:
                st["r2list"].append(d)
            rows.append(row)
            continue
        q = Decimal(16) * Vbsd / (fact * b["om1"] * reg * b["tam"])
        k = max(int(q.sqrt().to_integral_value(rounding=ROUND_FLOOR)), 1)
        # rank-2 L''(1) via finite differences carries ~1e-14 noise; others 40+ digits
        tol = Decimal("1e-11") if mode == "r2" else Decimal("1e-24")
        ok = abs(q - k * k) < tol * max(k * k, 1) or abs(q - (k + 1) ** 2) < tol * (k + 1) ** 2
        ks = k if abs(q - k * k) < abs(q - (k + 1) ** 2) else k + 1
        row.update({"mode": mode, "Reg": str(reg)[:30], "Q": f"{q:.6g}",
                    "sha_pred": ks if ok else "?", "Q_square": bool(ok)})
        if ok:
            st["sqc"] += 1
            dev = abs(q - ks * ks)
            md = st.setdefault("maxdev", {})
            if dev > md.get(mode, Decimal(0)):
                md[mode] = dev
            st["sha"][str(ks)] = st["sha"].get(str(ks), 0) + 1
            # engine-ball certified uniqueness of square (if ball decisive and r_an==rk)
            if decisive and r_an is not None and r_an == rk:
                qe_lo = Decimal(16) * (val - rad) / (fact * b["om1"] * reg * b["tam"])
                qe_hi = Decimal(16) * (val + rad) / (fact * b["om1"] * reg * b["tam"])
                klo = math.isqrt(int(qe_lo)) if qe_lo > 0 else 0
                khi = math.isqrt(int(qe_hi)) + 1
                sqs = [x * x for x in range(max(klo, 1), khi + 2) if qe_lo <= x * x <= qe_hi]
                if len(sqs) == 1:
                    st["certsq"] += 1
        else:
            st["nonsq"].append((d, mode, f"{q:.8g}"))
        rows.append(row)

    keys = ["d", "N", "w", "w_formula", "w_rootno", "pari_rank", "mode", "Reg",
            "Q", "sha_pred", "Q_square"]
    with open(os.path.join(HERE, "results.csv"), "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        wr.writeheader()
        for r in rows:
            wr.writerow(r)

    print("curves:", len(ds), "| engine certs:", len(certs), "| bsd:", len(bsd),
          "| rank:", len(rank), "| rootno:", len(rootno), "| reg2:", len(reg2))
    print("N=32d^2 failures:", st["Nbad"] or "none")
    print("root number disagreements (w_formula vs ellrootno vs engine):", st["wbad"] or "none")
    print("engine-vs-PARI value disagreements:", st["vdisagree"] or "none",
          "| maxdiff(<=1e-30 cap):", st["maxvd"])
    print("r_an certified: rank0 =", st["ran"][0], ", rank1 =", st["ran"][1],
          ", ball-not-decisive =", st["ran"][None])
    print("PARI rank vs r_an: agree", st["rankagree"], "| disagree:", st["rankdis"] or "none")
    print("rank-2 curves (rk=2):", len(st["r2list"]), "d =", st["r2list"])
    if st.get("nogen"):
        print("rank-1 without PARI generator (Reg unknown), n =", len(st["nogen"]), ":", st["nogen"][:50])
    if st["reg2miss"]:
        print("rank-2 missing reg data:", st["reg2miss"])
    print("Q(d) perfect squares:", st["sqc"], "(engine-ball certified uniquely:",
          st["certsq"], ") | non-square outliers:", st["nonsq"] or "none")
    print("predicted #Sha distribution:", dict(sorted(st["sha"].items(), key=lambda kv: int(kv[0]))))
    if st.get("maxdev"):
        print("max |Q - k^2| per mode:", {k: str(v) for k, v in st["maxdev"].items()})
    if st["norank"]:
        print("no PARI rank data for d:", st["norank"][:30])


if __name__ == "__main__":
    main()
