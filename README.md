# A Size Law for Periodic Orbits of the Fibonacci–gcd Cofactor Map, and All Orbits up to Period 40: Reproduction Code

This repository reproduces every computation in

> Hawarey, M. (2026). A Size Law for Periodic Orbits of the Fibonacci–gcd Cofactor Map, and All Orbits up to Period 40. *AIR Journal of Mathematics and Computational Sciences*.

## The objects

- **The recurrence.** C(n) = C(n−1) + C(n−2) + gcd(C(n−1), C(n−2)).
- **The cofactor map.** Its gcd dynamics is the map R(α, β) = (β/g, (α+β+1)/g) on coprime pairs, where g = gcd(β, α+1).
- **Amplifications.** A step with g > 1 is an amplification. A cycle of period q has k amplifications and cycle multiplier G, the product of the step multipliers over one period. We write t = G − L_q and τ = G − φ^q.
- **The golden coordinate.** Y = (2u − v) + (v − u)φ with u = α+1, v = β+1.

## Usage

```
python3 paper7_reproduce.py                    # full run, about 20 seconds
python3 paper7_reproduce.py lucas product      # selected parts
python3 paper7_reproduce.py search --A 10000000   # the Python search up to alpha = 10^7 (about 2 minutes)
```

Python 3.8 or later is required, with no third-party packages. All arithmetic is exact (integers, rationals, and exact arithmetic in Q(√5)); decimal values at 150 digits are used only for printed ratios. Every check prints `[ok]` or `[FAIL]`, and the run ends with `ALL CHECKS PASSED` or a list of the failed checks. `paper7_full_run.log` is the output of a full run.

## Parts and the paper's labels

| Part | Reproduces |
|---|---|
| `product` | Lemma 1, Theorem 2 (formula (4), exactly in Q(√5)) and Corollary 3 (bounds (5)) on the 24 distinct cycles of Section 4, Proposition 7 and period 44; Figure 1 ratios Y*τ/G from 0.50 to 2.67 |
| `bound` | Proposition 5: A(q) for q ≤ 40, with A(40) = 4.598·10^9 < 5·10^9 |
| `lucas` | Proposition 7: all 490,207 patterns with G = L_q, even q ≤ 40, every k; 19 cycles |
| `search` | Theorem 6 in pure Python for α ≤ A (default 10^6): the same cycles as Table 1 in that range |
| `logs` | Theorem 6 from the archived logs of the full C run: 21 blocks, 3,039,635,462 starts, the six cycles of Table 1 |
| `theorem8` | Theorem 8 and Corollary 9, Table 1: exactly 22 cycles with q ≤ 40; amplification counts {1, 3, 6, 10} |
| `larger` | Proposition 10 and Table 2: the cycles with t ≥ 1, their multipliers and G/t |

## The full search of Theorem 6 (C)

The search to α = 5·10^9 took about 9.6 core-hours. It uses 128-bit integers (GCC or Clang on a 64-bit machine):

```
gcc -O2 -o paper7_search_original paper7_search_original.c -lm
./paper7_search_original 1 1000000 5000 > paper7_g_1.txt 2> paper7_g_1.err
```

The run used 21 blocks of α: [1, 10^6], [10^6+1, 2.5·10^8], and then blocks of 2.5·10^8 up to 5·10^9. The blocks are independent and can run in parallel. For each α the program takes the golden-line start of Lemma 1 (u = α+1, v = β+1 with v − uφ in (−1/φ, 1/φ²)), skips it unless gcd(α, β) = 1, and iterates R. It stops when the first entry drops below α, when an entry exceeds 2^100, or after 5,000 steps.

Log format, in `search_logs/`:

- `paper7_g_<A0>.txt`: one line per cycle found, `period k alpha beta`, where (α, β) is the least state of the cycle.
- `paper7_g_<A0>.err`: `tested N`, the number of starts tested in the block.

`paper7_search.c` is the same search with one output line per start, `alpha beta steps amplifications reason` (reason 0 = returned to the start, a cycle; 1 = first entry dropped below α; 2 = an entry above 2^100; 3 = 5,000 steps). Add `--cycles-only` to print only the cycles. Its per-start output was compared with an independent Python implementation on 7.3 million starts and agreed on every start.

```
gcc -O2 -o paper7_search paper7_search.c -lm
./paper7_search 1 1000000 5000 --cycles-only
```

## Citation

See `CITATION.cff`. License: MIT (see `LICENSE`).
