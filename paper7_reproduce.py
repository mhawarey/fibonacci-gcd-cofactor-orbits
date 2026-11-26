#!/usr/bin/env python3
"""paper7_reproduce.py: reproduces every number quoted in
   Hawarey, M. (2026). A Size Law for Periodic Orbits of the Fibonacci-gcd Cofactor Map, and All Orbits up to Period 40.
Pure Python 3.8+, no third-party packages. Labels are those of the paper.
Parts (run all by default, or name some):  product  bound  lucas  search  logs  theorem8  larger
  product  : Lemma 1, Theorem 2 (exact in Q(sqrt5)), Corollary 3 on the 24 distinct cycles of Section 4, Proposition 7 and q = 44
  bound    : Proposition 5, A(q) for q <= 40
  lucas    : Proposition 7, all patterns with G = L_q, even q <= 40, every k (under a minute)
  search   : Theorem 6 in pure Python for alpha <= A (default A = 10^6; option --A N); the full 5*10^9 run uses
             paper7_search_original.c (see README)
  logs     : Theorem 6 from the archived logs of the full run (search_logs/): starts tested and the cycles found
  theorem8 : Theorem 8 and Corollary 9, Table 1
  larger   : Proposition 10, Table 2
"""
import sys, os, glob, time, itertools
from math import gcd, isqrt, log2
from fractions import Fraction as Fr
from decimal import Decimal, getcontext
getcontext().prec = 150
R5 = Decimal(5).sqrt(); PHI = (1 + R5) / 2; PSI = 1 - PHI
HERE = os.path.dirname(os.path.abspath(__file__))
FAILS = []
def check(name, cond):
    print(("  [ok] " if cond else "  [FAIL] ") + name)
    if not cond: FAILS.append(name)

# ---------- basics ----------
F, L = [0, 1], [2, 1]
def fib(n):
    while len(F) <= n + 2: F.append(F[-1] + F[-2]); L.append(L[-1] + L[-2])
    return F[n]
def luc(n): fib(n); return L[n]
def R(s):
    a, b = s; g = gcd(b, a + 1); return (b // g, (a + b + 1) // g), g
def golden(s):  # Y = (2u - v) + (v - u) phi as integer pair
    u, v = s[0] + 1, s[1] + 1; return (2 * u - v, v - u)
def num(Y, c): return Decimal(Y[0]) + Decimal(Y[1]) * c
def mul(x, y): a, b = x; c, d = y; return (a * c + b * d, a * d + b * c + b * d)
def phipow(n): return (fib(n - 1), fib(n)) if n >= 1 else (1, 0)
def state_of(Y): a, b = Y; return (a + b - 1, a + 2 * b - 1)
def cycle(s0):
    s, seq = s0, []
    while True:
        n, g = R(s); seq.append((s, n, g)); s = n
        if s == s0: return seq
        if len(seq) > 100000: raise ValueError("not periodic")
def describe(s0):
    seq = cycle(s0); q = len(seq); G = 1
    for _, _, g in seq: G *= g
    posts = [(n, g) for _, n, g in seq if g > 1]; k = len(posts)
    idx = [i for i, (_, _, g) in enumerate(seq) if g > 1]
    segs = [(idx[(j + 1) % k] - idx[j]) % q or q for j in range(k)]
    amin = min(s[0] for s, _, _ in seq); smin = [s for s, _, _ in seq if s[0] == amin][0]
    return dict(q=q, k=k, G=G, t=G - luc(q), posts=posts, smin=smin, amin=amin,
                gs=[g for _, g in posts], segs=segs, seq=seq)

# ---------- Lucas patterns (Proposition 7) ----------
def factor(n):  # trial division; enough for L_q, q <= 40 (L_40 < 2.3e8), and a proof of primality
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
        p += 1 if p == 2 else 2
    if n > 1: f[n] = f.get(n, 0) + 1
    return f
def compositions(n, k):
    if k == 1: yield (n,); return
    for a in range(1, n - k + 2):
        for rest in compositions(n - a, k - 1): yield (a,) + rest
def ordered_factorizations(n, k, divs):
    if k == 1:
        if n >= 2: yield (n,)
        return
    for d in divs:
        if d >= 2 and n % d == 0 and n // d >= 2 ** (k - 1):
            for rest in ordered_factorizations(n // d, k - 1, divs): yield (d,) + rest
def canonical(qs, gs):
    return min(tuple(zip(qs[i:] + qs[:i], gs[i:] + gs[:i])) for i in range(len(qs)))
def pattern_test(qs, gs):
    """(R5) for G = L_q: Y_0 by the unit formula; content 1 before each amplification, exact multipliers,
    least period q by integer iteration of R. Returns the start state or None."""
    q = sum(qs); Z = (0, 0); P = 1; tail = q
    for qj, gj in zip(qs, gs):
        tail -= qj; t_ = mul((fib(qj - 1) - 1, fib(qj)), phipow(tail)); Z = (Z[0] + P * t_[0], Z[1] + P * t_[1]); P *= gj
    Y = mul(phipow(q), Z); Y = (Y[0] + 1, Y[1]); s0 = state_of(Y)
    if s0[0] < 1 or s0[1] < 1: return None
    W = Y
    for qj, gj in zip(qs, gs):
        w = W
        for i in range(qj):
            if gcd(w[0] - 1, w[1]) > 1: return None
            w = mul(w, (0, 1))
        if gcd(w[0] - 1, w[1]) != gj: return None
        W = (1 + (w[0] - 1) // gj, w[1] // gj)
    s = s0; amps = 0
    for j in range(q):
        s, m = R(s); amps += m > 1
        if s == s0 and j < q - 1: return None
    return s0 if (s == s0 and amps == len(qs)) else None
def lucas_census(QMAX=40, verbose=True):
    total, found = 0, []
    for q in range(2, QMAX + 1, 2):
        Lq = luc(q); f = factor(Lq); divs = [1]
        for p, e in f.items(): divs = [d * p ** i for d in divs for i in range(e + 1)]
        divs.sort(); Om = sum(f.values()); nq = 0
        for K in range(1, min(Om, q) + 1):
            seen = set()
            for gs in ordered_factorizations(Lq, K, divs):
                for qs in compositions(q, K):
                    c = canonical(qs, gs)
                    if c in seen: continue
                    seen.add(c); nq += 1; s0 = pattern_test(list(qs), list(gs))
                    if s0: found.append((q, K, qs, gs, s0))
        total += nq
        if verbose: print(f"    q = {q:2d}: L_q = {Lq} = {' * '.join(f'{p}^{e}' if e > 1 else str(p) for p, e in sorted(f.items()))}, patterns {nq}")
    return total, found

KNOWN_SEARCH = [(1, 2), (3, 5), (29, 47), (97, 157), (3870813, 6263107), (4003, 6477)]
P6_Q44 = (675056745204946765211, 1092264758076481467317)

def part_product(lucas_states):
    print("\n== Lemma 1, Theorem 2, Corollary 3 (Section 3)")
    cyc = KNOWN_SEARCH + lucas_states + [P6_Q44]; ratios = []
    for s0 in cyc:
        d = describe(s0); q, k, G = d["q"], d["k"], d["G"]; tau = Decimal(G) - PHI ** q
        ok = tau > 0
        for s, n, g in d["seq"]:
            Yp = num(golden(s), PSI); ok &= abs(Yp) < 1
            if g > 1: ok &= 1 - PHI / g < num(golden(n), PSI) < 1 - 1 / (PHI * PHI * g)
        ok &= any(n[0] == d["amin"] for n, g in d["posts"])
        u = d["amin"] + 1; ok &= d["smin"][1] + 1 == (u + isqrt(5 * u * u)) // 2
        # Theorem 2, exactly in Q(sqrt5): elements a + b*phi with rational a, b
        prod = (Fr(G), Fr(0))
        for n, g in d["posts"]:
            Y = golden(n); nrm = Y[0] * Y[0] + Y[0] * Y[1] - Y[1] * Y[1]; inv = (Fr(Y[0] + Y[1], nrm), Fr(-Y[1], nrm))
            x = (Fr(g - 1, g) * inv[0], Fr(g - 1, g) * inv[1]); prod = mul(prod, (1 - x[0], -x[1]))
        ok &= prod == (Fr(fib(q - 1)), Fr(fib(q)))
        pc = Decimal(G)
        for n, g in d["posts"]: pc *= 1 - (1 - Decimal(1) / g) / num(golden(n), PSI)
        ok &= abs(pc - PSI ** q) < Decimal(10) ** -80
        Ys = min((num(golden(n), PHI), g) for n, g in d["posts"])
        ok &= (1 - Decimal(1) / Ys[1]) * G / tau <= Ys[0] < k * G / tau
        ok &= Decimal("0.3618") * G / tau < d["amin"] + 1 < (PHI / R5) * k * G / tau + Decimal("0.2764")
        ratios.append(float(Ys[0] * tau / G))
        check(f"q = {q:2d}, k = {k:2d}, t = {d['t']}: Lemma 1, (4) exact, (5); Y*tau/G = {float(Ys[0] * tau / G):.3f}", ok)
    print(f"  cycles checked: {len(cyc)}; Y*tau/G between {min(ratios):.2f} and {max(ratios):.2f}")
    print(f"  k = 1: Y*tau/G = 1/2, 2/3, 6/7 -> {[round(r, 4) for r in ratios[:3]]}")
    print(f"  constants: phi/sqrt5 = {float(PHI / R5):.5f} (< 0.7236), -psi/sqrt5 = {float(-PSI / R5):.5f} (< 0.2764), |psi| = {float(-PSI):.4f}")

def part_bound():
    print("\n== Proposition 5")
    A = lambda q: (PHI / R5) * Decimal(log2(luc(q) + 1)) * (luc(q) + 1) / (1 + PSI ** q) + Decimal("0.2764")
    vals = [A(q) for q in range(1, 41)]
    check("A(q) increases for q <= 40", all(vals[i] < vals[i + 1] for i in range(39)))
    for q in (36, 37, 40): print(f"  A({q}) = {float(A(q)):.3e}")
    check("A(40) < 5*10^9", A(40) < 5 * 10 ** 9)
    check("second branch: q + 2.1 below the first, q = 1..40", all(Decimal(q) + Decimal("2.1") < (A(q) - Decimal("0.2764")) / (PHI / R5) for q in range(1, 41)))
    Ystar = Decimal(5 * 10 ** 9 + 1) / (PHI / R5); top = Decimal("1.18") * PHI ** 40 * Ystar
    print(f"  Y* < {float(Ystar):.3e} = {float(Ystar / (5 * 10 ** 9)):.3f} * 5e9; entries < {float(top):.2e} < 2e18 < 2^100")
    check("entries below 2*10^18 for q <= 40", top < 2 * 10 ** 18)

def part_lucas():
    print("\n== Proposition 7 (G = L_q, even q <= 40, every k)")
    t0 = time.time(); total, found = lucas_census(40)
    print(f"  patterns: {total:,}; cycles: {len(found)} [{time.time() - t0:.0f} s]")
    from collections import Counter
    print("  by (q, k):", sorted(Counter((f[0], f[1]) for f in found).items()))
    a28 = sorted(describe(f[4])["amin"] for f in found if f[0] == 28); a34 = [describe(f[4])["amin"] for f in found if f[0] == 34]
    print(f"  least first entries: q = 28 from {a28[0]:.2e} to {a28[-1]:.2e}; q = 34: {a34[0]:.2e}")
    check("490,207 patterns and 19 cycles", total == 490207 and len(found) == 19)
    return [f[4] for f in found if f[1] == 3]

def part_search(A):
    print(f"\n== Theorem 6 in Python for alpha <= {A:,} (the full 5*10^9 run: paper7_search_original.c, see README)")
    t0 = time.time(); n = 0; cyc = []; M = 1 << 100
    for a in range(1, A + 1):
        u = a + 1; v = (u + isqrt(5 * u * u)) // 2
        if not (5 * (u - 1) ** 2 < (2 * (v - 1) - (u - 1)) ** 2): v += 1
        b = v - 1
        if b < 1 or gcd(a, b) != 1: continue
        n += 1; s = (a, b)
        for i in range(1, 5001):
            s, g = R(s)
            if s[1] > M or s[0] < a: break
            if s == (a, b): cyc.append((a, b)); break
    print(f"  starts {n:,}; cycles through their least state: {cyc} [{time.time() - t0:.0f} s]")
    check("same cycles as Table 1 within this range", cyc == [c for c in sorted(KNOWN_SEARCH) if c[0] <= A])

def part_logs():
    print("\n== Theorem 6 from the archived logs of the full run (search_logs/)")
    errs = sorted(glob.glob(os.path.join(HERE, "search_logs", "*.err"))); outs = sorted(glob.glob(os.path.join(HERE, "search_logs", "*.txt")))
    starts = sum(int(open(f).read().split()[1]) for f in errs)
    ranges = sorted(int(os.path.basename(f).split("_")[-1].split(".")[0]) for f in errs)
    found = set()
    for f in outs:
        for line in open(f):
            p = line.split()
            if len(p) == 4: found.add((int(p[2]), int(p[3])))
    print(f"  blocks: {len(errs)} starting at {ranges}; starts tested: {starts:,}")
    check("21 blocks covering alpha = 1 .. 5*10^9", len(errs) == 21 and ranges[0] == 1 and ranges[-1] == 4750000001)
    check("starts tested 3,039,635,462", starts == 3039635462)
    check("cycles found = the six of Table 1", found == set(KNOWN_SEARCH))
    print(f"  cycles per start: {6 / starts:.2e} (about two per billion)")

def part_theorem8(lucas_states):
    print("\n== Theorem 8, Corollary 9, Table 1")
    rows = sorted((describe(s) for s in KNOWN_SEARCH), key=lambda d: d["q"])
    print("  Table 1: q, k, t, state, Y*tau/G")
    for d in rows:
        Ys = min((num(golden(n), PHI), g) for n, g in d["posts"]); tau = Decimal(d["G"]) - PHI ** d["q"]
        print(f"   {d['q']:3d} {d['k']:3d} {d['t']:>15,} {d['smin']} {float(Ys[0] * tau / d['G']):.3f}")
    small = [d for d in rows if d["q"] <= 40 and d["t"] >= 1]
    luc_ = [describe(s) for s in lucas_states] + [describe((3, 5)), describe((29, 47))]
    allc = small + luc_
    from collections import Counter
    byk = Counter(d["k"] for d in allc); least = {}
    for d in allc: least[d["k"]] = min(least.get(d["k"], 999), d["q"])
    print(f"  cycles with q <= 40: {len(allc)}; by k: {sorted(byk.items())}; least period by k: {sorted(least.items())}")
    check("Theorem 8: 22 cycles", len(allc) == 22)
    check("Corollary 9: k in {1, 3, 6, 10}; least periods 28 (k = 3), 22 (k = 6), 39 (k = 10)",
          set(byk) == {1, 3, 6, 10} and least[3] == 28 and least[6] == 22 and least[10] == 39)

def part_larger():
    print("\n== Proposition 10, Table 2")
    for s in KNOWN_SEARCH:
        d = describe(s)
        if d["t"] < 1: continue
        f = factor(d["G"]) if d["G"] < 10 ** 16 else None
        fs = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(f.items())) if f else ""
        print(f"  q = {d['q']}, k = {d['k']}, G = {d['G']:,} = {fs}, t = {d['t']:,}, alpha_min = {d['amin']:,}, G/t = {d['G'] / d['t']:.4g}")
        if d["k"] >= 2:
            print(f"     g = {tuple(d['gs'])}, q_j = {tuple(d['segs'])}; threes: {d['gs'].count(3)} of {d['k']}; q > 2k: {d['q'] > 2 * d['k']}; 15 | G: {d['G'] % 15 == 0}")
            check(f"q = {d['q']}: all multipliers odd and q > 2k", d["q"] > 2 * d["k"] and all(g % 2 for g in d["gs"]))

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    A = int(sys.argv[sys.argv.index("--A") + 1]) if "--A" in sys.argv else 10 ** 6
    if "--A" in sys.argv: args = [a for a in args if a != str(A)]
    parts = args or ["lucas", "product", "bound", "search", "logs", "theorem8", "larger"]
    t0 = time.time(); lucas_states = None
    if "lucas" in parts or "product" in parts or "theorem8" in parts:
        lucas_states = part_lucas() if "lucas" in parts else [f[4] for f in lucas_census(40, verbose=False)[1] if f[1] == 3]
    if "product" in parts: part_product(lucas_states)
    if "bound" in parts: part_bound()
    if "search" in parts: part_search(A)
    if "logs" in parts: part_logs()
    if "theorem8" in parts: part_theorem8(lucas_states)
    if "larger" in parts: part_larger()
    print(f"\nTotal time {time.time() - t0:.0f} s.", "ALL CHECKS PASSED" if not FAILS else f"FAILED: {FAILS}")
