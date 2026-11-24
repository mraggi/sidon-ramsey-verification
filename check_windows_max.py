#!/usr/bin/env python3
"""Check that every window in windows.json is the longest one, using only Python's standard library.

    python3 check_windows_max.py

verify_windows.py checks that each stored (u, t, m) gives a Sidon k-partition of [v+m].
This script checks the two other columns of the table of windows (Appendix B):

  * m is the longest window over all affine images x -> u*x + t (u a unit of Z_v);
  * m_0 is the largest m with k*m*(m-1) < 2v (Proposition 4.3).

Reading the image x -> u*x + t from s = 1 colours x by the class of u*x + t, so the window
is the run of residues t+u, t+2u, ..., t+mu. For a fixed u, changing t only shifts the
cyclic sequence class(u*i), i in Z_v, so the longest window for u is the longest cyclic
run of that sequence whose residues lie in different classes. For Singer partitions the
residue 0 lies on every line and may go in any class not used by the rest of the window,
so it never conflicts; a window still has at most k residues.
"""
import json
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent
WILD = -1  # the residue 0 of a Singer partition


def labels(e):
    v = e["v"]
    if e["construction"] == "Singer":
        line = set(e["line"])
        lab = [WILD] * v
        for x in range(1, v):
            ds = [d for d in line if (x + d) % v in line]
            assert len(ds) == 1, "not a perfect difference set"
            lab[x] = ds[0]
        return lab
    p, g = e["p"], e["g"]
    pw = [pow(g, j, p) for j in range(p - 1)]
    return [(n % p - pw[n % (p - 1)]) % p for n in range(v)]


def longest_window(e):
    v, k = e["v"], e["k"]
    lab = labels(e)
    best = 0
    for u in range(1, v):
        if gcd(u, v) != 1:
            continue
        seq = [lab[(u * i) % v] for i in range(v)]
        seq += seq[:k]
        count = {}
        lo = 0
        for hi, x in enumerate(seq):
            if x != WILD:
                count[x] = count.get(x, 0) + 1
                while count[x] > 1:
                    y = seq[lo]
                    lo += 1
                    if y != WILD:
                        count[y] -= 1
            best = max(best, min(hi - lo + 1, k))
    return best


def main():
    rows = json.loads((HERE / "windows.json").read_text())
    bad = 0
    for e in rows:
        k, v, m = e["k"], e["v"], e["m"]
        m0 = max(n for n in range(1, v) if k * n * (n - 1) < 2 * v)
        mx = longest_window(e)
        ok = mx == m
        bad += not ok
        print(f"{'OK  ' if ok else 'FAIL'} {e['construction']:6} k={k:2} v={v:4} "
              f"m={m:2} longest={mx:2} m0={m0:2}")
    print(f"{len(rows) - bad} of {len(rows)} windows are the longest")
    raise SystemExit(1 if bad else 0)


if __name__ == "__main__":
    main()
