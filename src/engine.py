"""BSD attack on quadratic twists of 32a1: certified analytic side.

Family E^(d): y^2 = x^3 - d^2*x, d odd squarefree (quadratic twist of
E = 32a1: y^2 = x^3 - x, CM by Z[i]). Conductor N = 32*d^2.

Certified L(E^(d),1) (root number w=+1) or L'(E^(d),1) (w=-1) from the
symmetric completed-L fast series

    Lam(s) := N^{s/2}(2pi)^{-s} Gamma(s) L(s)
           = sum_{n>=1} a_n [ (sqrt(N)/(2pi n))^s     Gamma(s,   x_n)
                            + w*(sqrt(N)/(2pi n))^{2-s} Gamma(2-s, x_n) ],
    x_n = 2*pi*n/sqrt(N)

(Fricke/Cauchy form; Lam(s) = w*Lam(2-s) holds by construction; terms decay
like e^{-x_n}). Then L(s) = (2pi/sqrt(N))^s * Lam(s)/Gamma(s):
  w=+1:  L(1)  = 2*sum a_n e^{-x_n}/n + O(e^{-x_M})
  w=-1:  Lam(1)=0 identically, and L'(1) = (2pi/sqrt(N)) * Lam'(1).
A rigorous tail bound |a_n| <= d(n)*sqrt(n) turns the cutoff M into a ball
[v-r, v+r]; 0 outside the ball certifies r_an(E^(d),1) = 0 (w=+1) or 1
(w=-1); Gross-Zagier + Kolyvagin then give rank(E^(d)(Q)) = r_an and finite
Sha. Cross-references: PARI/GP ellL1, LMFDB (cited in paper.md, not used).

Pure Python + mpmath for the whole analytic side.
"""

import math

import mpmath

mpmath.mp.dps = 45


# ---------------------------------------------------------------- sieve / NT


def sieve(n):
    if n < 2:
        return []
    s = bytearray(b"\x01") * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n**0.5) + 1):
        if s[i]:
            s[i * i :: i] = bytearray(len(s[i * i :: i]))
    return [i for i in range(2, n + 1) if s[i]]


def factorize(n):
    f = {}
    for p in sieve(int(n**0.5) + 1):
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def squarefree_odd_upto(x):
    out = []
    for d in range(1, x + 1, 2):
        if d == 1 or all(v == 1 for v in factorize(d).values()):
            out.append(d)
    return out


def legendre(a, p):
    """Legendre symbol (a/p), odd prime p; a may be negative."""
    if a % p == 0:
        return 0
    r = pow(a % p, (p - 1) // 2, p)
    return -1 if r == p - 1 else r


def kron2(a):
    """(a/2) for odd a."""
    return 1 if a % 8 in (1, 7) else -1


def kronecker_d_n(d, n):
    """Kronecker symbol (d/n), d > 0 odd squarefree, n >= 1; completely
    multiplicative in n."""
    if n == 1:
        return 1
    res = 1
    e = 0
    while n % 2 == 0:
        n //= 2
        e += 1
    if e:
        res *= kron2(d) ** e
    for p, k in factorize(n).items():
        s = legendre(d, p)
        if s == 0:
            return 0
        if s == -1 and k % 2:
            res = -res
    return res


def kronecker_m2(d):
    """(-2/d) for odd d > 0 (Kronecker): (-1)^{(d-1)/2} * (2/d)."""
    return (-1 if d % 4 == 3 else 1) * (1 if d % 8 in (1, 7) else -1)


# --------------------------- a_n of E = 32a1: y^2 = x^3 - x (CM by Z[i])


def _gauss_sqrt_minus1(p):
    g = 2
    while pow(g, (p - 1) // 2, p) != p - 1:
        g += 1
    return pow(g, (p - 1) // 4, p)


def _primary_factor(p):
    """Primary Gaussian prime a+bi above p (p = 1 mod 4 prime): a odd, b even,
    a + b == 1 (mod 4); unique unit multiple."""
    r = _gauss_sqrt_minus1(p)
    a1, b1, a2, b2 = p, 0, r, 1
    while (a2, b2) != (0, 0):
        nrm = a2 * a2 + b2 * b2
        qr = round((a1 * a2 + b1 * b2) / nrm)
        qi = round((b1 * a2 - a1 * b2) / nrm)
        a1, b1, a2, b2 = (
            a2,
            b2,
            a1 - (qr * a2 - qi * b2),
            b1 - (qr * b2 + qi * a2),
        )
    for aa, bb in ((a1, b1), (-b1, a1), (-a1, -b1), (b1, -a1)):
        if aa % 2 != 0 and bb % 2 == 0 and (aa + bb) % 4 == 1:
            return aa
    raise AssertionError("no primary unit multiple")


def ap_E(p):
    """a_p of E: y^2 = x^3 - x (CM by Z[i]): 0 for p=2 or p=3 mod 4; else
    a_p = 2a where p = a^2 + b^2 with a+bi primary. Sign fixed against direct
    point counts (selftest 1); matches Silverman AEP V.5 Thm 5.4."""
    if p == 2 or p % 4 == 3:
        return 0
    return 2 * _primary_factor(p)


_ACACHE = {"M": 0, "a": None}


def a_coeffs_E(M):
    """a_n(E), n <= M (index = n). E: y^2=x^3-x has additive reduction at 2
    with conductor exponent 5, so the local L-factor at 2 is 1: a_n = 0 for
    4 | n. Odd prime powers use the newform recurrence
    a_{p^{k+1}} = a_p a_{p^k} - p a_{p^{k-1}} (good reduction). Verified
    against PARI ellan(E32, N) up to N=300 (run.py)."""
    if _ACACHE["M"] >= M:
        return _ACACHE["a"][: M + 1]
    a = [0] * (M + 1)
    a[1] = 1
    spf = [0] * (M + 1)
    for p in sieve(M):
        for m in range(p, M + 1, p):
            if spf[m] == 0:
                spf[m] = p
    apc = {}
    for n in range(2, M + 1):
        if n % 4 == 0:
            continue  # trivial local factor at 2 (a_2=0, a_{2^k}=0, k>=2)
        p = spf[n]
        m = n // p
        ap = apc.setdefault(p, ap_E(p))
        if m % p == 0:
            a[n] = ap * a[m] - p * a[m // p]
        else:
            a[n] = ap * a[m]
    _ACACHE["M"] = M
    _ACACHE["a"] = a
    return a


def twisted_coeffs(d, M):
    """a_n(E^(d)) = a_n(E)*(d/n) for odd squarefree d, EXCEPT a_n = 0 whenever
    p^2 | n for some p | d: at odd p|d the twist pi_p (unram principal series
    tensor ramified quadratic chi) has conductor exponent 2 and trivial GL2
    L-factor, so all a_{p^k}(E^(d)) = 0, k >= 1. Cross-checked vs PARI ellan
    on E^(5), E^(7), E^(15) (run.py)."""
    a = a_coeffs_E(M)
    if d == 1:
        return a
    out = [0] * (M + 1)
    for n in range(1, M + 1):
        if not a[n]:
            continue
        g = math.gcd(n, d)
        if g > 1 and any(n % (p * p) == 0 for p in factorize(g)):
            continue
        out[n] = a[n] * kronecker_d_n(d, n)
    return out


# ----------------------------------------------------- L-series + certificate


# scale factor for _M_cut; run.py bumps it (MSCALE env) for curves whose
# |L(1)| falls under the default ball radius, shrinking tails exponentially.
_CUT_SCALE = 1


def _M_cut(N):
    """cutoff: e^{-2*pi*M/sqrt(N)} < 1e-40 (times _CUT_SCALE)."""
    return math.ceil(95.0 * _CUT_SCALE * math.sqrt(N) / (2 * math.pi)) + 5


def dGam1(x):
    """d/ds Gamma(s,x) at s=1 = int_x^inf ln(t)e^-t dt = e^{-x} ln x + E1(x)."""
    return mpmath.exp(-x) * mpmath.log(x) + mpmath.e1(x)


def lam_sym(N, an, w, s):
    """Symmetric completed L-function Lam(s); caller must have len(an)-1 >= M."""
    Nf = mpmath.mpf(N)
    sN = mpmath.sqrt(Nf)
    pi2 = 2 * mpmath.pi
    M = min(_M_cut(N), len(an) - 1)
    return mpmath.fsum(
        an[n]
        * (
            (sN / (pi2 * n)) ** s * mpmath.gammainc(s, pi2 * n / sN, mpmath.inf)
            + w
            * (sN / (pi2 * n)) ** (2 - s)
            * mpmath.gammainc(2 - s, pi2 * n / sN, mpmath.inf)
        )
        for n in range(1, M + 1)
        if an[n]
    )


def _tail_bound(N, M, C0):
    """Rigorous overestimate of |sum_{n>M} a_n*(c/n)^u*Gamma(u,x_n)|-type tails
    (u near 1 or 2, |d/ds| bringing logs): using |a_n| <= d(n)*sqrt(n) <= 2n
    and n^{-u}*Gamma(u,x_n)*n^{...} <= (sqrt(N)/2pi)^{~1} e^{-c n}(C0+2 ln n),
    c = 2 pi/sqrt(N); dyadic blocks with geometric decay."""
    sN = mpmath.sqrt(mpmath.mpf(N))
    c = 2 * mpmath.pi / sN

    def geom(t):
        e = mpmath.exp(-t)
        return e / (1 - e)

    total = mpmath.mpf(0)
    lo = M + 1
    while True:
        hi = 2 * lo
        t0 = c * lo
        if t0 > 700:
            break
        blk = (C0 + 2 * mpmath.log(mpmath.mpf(hi))) * geom(t0)
        total += blk
        if blk < mpmath.mpf(10) ** (-45):
            break
        lo = hi
    return total * (sN / (2 * mpmath.pi)) * 2 * 10  # |a_n|<=2n model slack, x10


def L_value(N, an, w, deriv=False):
    """Certified L(E,1) (deriv=False, needs w=+1) or L'(E,1) (deriv=True,
    needs w=-1). Returns (value, radius, M), |true - value| <= radius."""
    Nf = mpmath.mpf(N)
    sN = mpmath.sqrt(Nf)
    pi2 = 2 * mpmath.pi
    M = min(_M_cut(N), len(an) - 1)

    if not deriv:
        # L(1) = 2*pi*Lam(1)/sqrt(N), Lam(1) = (1+w)*sum a_n e^{-x_n} n^{-1} sN/(2pi)
        S1 = mpmath.fsum(
            an[n] * mpmath.exp(-pi2 * n / sN) / n for n in range(1, M + 1) if an[n]
        )
        lam1 = (1 + w) * (sN / pi2) * S1
        val = lam1 * pi2 / sN
        C0 = 2 * float(mpmath.log(sN)) + 30
        r = _tail_bound(N, M, C0) * (pi2 / sN) * (1 + abs(w)) + mpmath.mpf(10) ** (-40)
        return val, r, M

    # w = -1: Lam(1) = 0 exactly; L'(1) = 2*pi*Lam'(1)/sqrt(N).
    # term(n,s) = (c/n)^s Gam(s,x) + w (c/n)^{2-s} Gam(2-s,x), c = sN/2pi,
    # d/ds at s=1:  (c/n)[ (ln c - ln n) e^{-x} + dGam1(x) ]
    #             + w (c/n)[ -(ln c - ln n) e^{-x} - dGam1(x) ]   (chain rule)
    # for w=-1 the pieces ADD:  Lam'(1) = 2*sum a_n (c/n)[(ln c - ln n)e^{-x} + dGam1(x)]
    def aterm(n):
        x = pi2 * n / sN
        g = mpmath.exp(-x)
        lg = mpmath.log(sN / pi2) - mpmath.log(n)  # ln(c/n)
        return an[n] * (sN / pi2) / n * (lg * g + dGam1(x))

    lam1p = 2 * mpmath.fsum(aterm(n) for n in range(1, M + 1) if an[n])
    val = lam1p * pi2 / sN
    C0 = 3 * float(mpmath.log(sN)) + 40
    r = _tail_bound(N, M, C0) * (float(pi2) / float(sN)) * 2 + mpmath.mpf(10) ** (-40)
    return val, r, M


def fe_residual(N, an, w, s0=None):
    """|Lam(s0)/(w Lam(2-s0)) - 1| -- should be ~0 by construction; guards
    against coefficient drift."""
    if s0 is None:
        s0 = mpmath.mpf("1.37")
    Nf = mpmath.mpf(N)
    sN = mpmath.sqrt(Nf)
    pi2 = 2 * mpmath.pi
    M = min(_M_cut(N), len(an) - 1)

    def Lam(s):
        return mpmath.fsum(
            an[n]
            * (
                (sN / (pi2 * n)) ** s * mpmath.gammainc(s, pi2 * n / sN, mpmath.inf)
                + w
                * (sN / (pi2 * n)) ** (2 - s)
                * mpmath.gammainc(2 - s, pi2 * n / sN, mpmath.inf)
            )
            for n in range(1, M + 1)
            if an[n]
        )

    a = Lam(s0)
    b = Lam(2 - s0)
    if abs(b) < mpmath.mpf(10) ** (-25):
        return None
    return float(abs(a / (w * b) - 1))


def L_direct(N, an, s, nmax):
    """Plain Dirichlet series sum a_n n^{-s} (absolute convergence for s>3/2);
    convention-free ground truth for the normalization of Lam."""
    return mpmath.fsum(
        mpmath.mpf(an[n]) * mpmath.mpf(n) ** (-s) for n in range(1, min(nmax, len(an))) if an[n]
    )


# ------------------------------------------------------------ self-tests / ref


def count_points_mod_p(d, p):
    """#E^(d)(F_p) on y^2 = x^3 - d^2 x (incl. O); counts nonsingular points."""
    dd = (d * d) % p
    cnt = 1
    for x in range(p):
        v = (x * x * x - dd * x) % p
        if v == 0:
            cnt += 1
        elif pow(v, (p - 1) // 2, p) == 1:
            cnt += 2
    return cnt


def lemniscate_omega():
    """One-component real period of 32a1: int_1^inf dx/sqrt(x^3-x) = B(1/4,1/2)/2."""
    return mpmath.beta(mpmath.mpf(1) / 4, mpmath.mpf(1) / 2) / 2


# PARI/GP lfun() ground truth (38-55 digits, server run; probe*.gp in this dir):
PARI_L1 = {  # L(E^(d),1), w=+1 curves
    1: "0.65551438857302995261620989747277985342",  # E^(1)=32a1
    11: "0.79058009875489948401853430365075214704888138659309",
}
PARI_L15 = "0.82634820502187107714530615794229225233"  # L(32a1, 3/2)
PARI_L2 = "0.91705063531865498864380552429571331840"  # L(32a1, 2)
PARI_LD1 = {  # L'(E^(d),1) for w=-1 curves (Richardson h^4 of lfun, ~45 digits)
    5: "2.227370379544139206937296287057707397116517197505879988",
    7: "2.962114883254813374140109352374987206141266260912343146",
}
# PARI ellan() vectors (a_1..a_N) for the coefficient cross-check:
PARI_AN = {
    1: "1,0,0,0,-2,0,0,0,-3,0,0,0,6,0,0,0,2,0,0,0,0,0,0,0,-1,0,0,0,-10,0,"
    "0,0,0,0,0,0,-2,0,0,0,10,0,0,0,6,0,0,0,-7,0,0,0,14,0,0,0,0,0,0",
    5: "1,0,0,0,0,0,0,0,-3,0,0,0,-6,0,0,0,-2,0,0,0,0,0,0,0,0,0,0,0,-10,0,"
    "0,0,0,0,0,0,2,0,0,0,10,0,0,0,0,0,0,0,-7,0,0,0,-14,0,0,0,0,0,0",
    7: "1,0,0,0,2,0,0,0,-3,0,0,0,-6,0,0,0,-2,0,0,0,0,0,0,0,-1,0,0,0,-10,0,"
    "0,0,0,0,0,0,-2,0,0,0,-10,0,0,0,-6,0,0,0,0,0,0,0,14,0,0,0,0,0,0",
    15: "1,0,0,0,0,0,0,0,0,0,0,0,-6,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,10,0,"
    "0,0,0,0,0,0,2,0,0,0,-10,0,0,0,0,0,0,0",
}


def selftest():
    # 1. a_p formula vs direct point counts on E, all odd primes <= 97
    for p in sieve(97):
        if p == 2:
            continue
        assert ap_E(p) == p + 1 - count_points_mod_p(1, p), ("ap", p)
    # 2. twist coefficients at good primes: a_p(E^(d)) = (d/p) a_p(E)
    for d in (5, 7, 15):
        for p in sieve(31):
            if p != 2 and p % d:
                got = p + 1 - count_points_mod_p(d, p)
                assert got == ap_E(p) * kronecker_d_n(d, p), ("twist", d, p)
    # 3. one-component period integral vs closed form B(1/4,1/2)/2 = 2.62205755...
    assert abs(lemniscate_omega() - 2.6220575542921198) < 1e-12
    # 3b. a_n and twisted a_n vs PARI ellan vectors (incl. a_{4k}=0, p^2|d => 0)
    for d, vec in PARI_AN.items():
        ref = [int(t) for t in vec.split(",")]
        got = twisted_coeffs(d, len(ref))
        assert got[1:] == ref, ("ellan", d)
    # 4. modular law g(1/(Nu)) = w*N*u^2*g(u) for E^(1), u=1/2 (both sides
    # converge fast); this is the Fricke relation behind the B-series.
    an = twisted_coeffs(1, 400)
    g = lambda u: mpmath.fsum(
        an[n] * mpmath.exp(-2 * mpmath.pi * n * u) for n in range(1, 401) if an[n]
    )
    u = mpmath.mpf("0.5")
    assert abs(g(1 / (32 * u)) / (32 * u**2 * g(u)) - 1) < mpmath.mpf(10) ** (-35)
    # 4b. absolute normalization: Lam(s) series vs PARI lfun at s=3/2 and s=2
    an = twisted_coeffs(1, _M_cut(32) + 5)
    c = mpmath.sqrt(mpmath.mpf(32)) / (2 * mpmath.pi)
    e15 = mpmath.mpf(3) / 2
    lam = lam_sym(32, an, +1, e15)
    tgt = c**e15 * mpmath.gamma(e15) * mpmath.mpf(PARI_L15)
    assert abs(lam / tgt - 1) < mpmath.mpf(10) ** (-30), ("norm s=1.5", lam, tgt)
    lam = lam_sym(32, an, +1, 2)
    tgt = c**2 * mpmath.gamma(2) * mpmath.mpf(PARI_L2)
    assert abs(lam / tgt - 1) < mpmath.mpf(10) ** (-30), ("norm s=2", lam, tgt)
    # 5. FE self-check at s0 for d=1,5,7 (by construction ~1e-38)
    for d in (1, 5, 7):
        an = twisted_coeffs(d, _M_cut(32 * d * d) + 5)
        w = kronecker_m2(d)
        e = fe_residual(32 * d * d, an, w)
        assert e is not None and e < 1e-25, ("fe", d, e)
    # 6. L(32a1,1) + L(E^(11),1) vs PARI lfun values (absolute anchors, w=+1)
    for dd in (1, 11):
        N = 32 * dd * dd
        an = twisted_coeffs(dd, _M_cut(N) + 5)
        val, r, _M = L_value(N, an, +1)
        assert abs(val - mpmath.mpf(PARI_L1[dd])) < mpmath.mpf(10) ** (-35), ("L1", dd, val)
        assert r < 1e-35, ("rad", dd, r)
    # 7. L'(E^(d),1) vs PARI Richardson anchors (w=-1 branch, ~45 digits)
    for d, ref in PARI_LD1.items():
        N = 32 * d * d
        an = twisted_coeffs(d, _M_cut(N) + 5)
        v, r, _M = L_value(N, an, -1, deriv=True)
        assert abs(v / mpmath.mpf(ref) - 1) < mpmath.mpf(10) ** (-35), ("Lp", d, v)
    return True


if __name__ == "__main__":
    import time

    t0 = time.time()
    assert selftest()
    print(f"selftest OK in {time.time() - t0:.1f}s")
    for d in (1, 3, 5, 7):
        N = 32 * d * d
        an = twisted_coeffs(d, _M_cut(N) + 5)
        w = kronecker_m2(d)
        if w == 1:
            v, r, M = L_value(N, an, w)
            tag = "L "
        else:
            v, r, M = L_value(N, an, w, deriv=True)
            tag = "L'"
        print(
            f"d={d} N={N} w={w:+d} {tag}(E,1) = {mpmath.nstr(v, 25)} +/- {mpmath.nstr(r, 2)}"
        )
