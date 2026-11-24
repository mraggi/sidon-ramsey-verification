#!/usr/bin/env python3
"""Check the fractional tiling of [108] by 12-element Sidon sets (paper, end of Section 3).

    python3 check_tiling108.py

tiling108.json lists 108 twelve-element subsets S of [108] with nonnegative rational
weights y_S. The script checks, in exact rational arithmetic, that every S is a Sidon set
and that sum_{S containing x} y_S = 1 for every x in [108], so that sum_S y_S = 9.
Consequently, if integer weights w give every 12-element Sidon subset of [108] weight at
least D, then sum_x w(x) = sum_S y_S w(S) >= 9D: the weight method of the paper cannot
exclude nine pairwise disjoint 12-element Sidon sets in [108].
"""
import json
from fractions import Fraction
from pathlib import Path

d = json.loads((Path(__file__).resolve().parent / "tiling108.json").read_text())
n, m = d["n"], d["m"]
sets, weights = d["sets"], [Fraction(w) for w in d["weights"]]
assert (n, m) == (108, 12) and len(sets) == len(weights)
cover = [Fraction(0)] * (n + 1)
for s, y in zip(sets, weights):
    assert y >= 0 and len(set(s)) == m and all(1 <= x <= n for x in s)
    sums = [a + b for i, a in enumerate(s) for b in s[i:]]
    assert len(set(sums)) == len(sums), f"not a Sidon set: {s}"
    for x in s:
        cover[x] += y
assert all(cover[x] == 1 for x in range(1, n + 1)), "not a tiling"
assert sum(weights) == 9
print(f"OK: {len(sets)} Sidon sets of size {m} with weights covering every point of [{n}] exactly once; total weight 9")
