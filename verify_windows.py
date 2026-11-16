#!/usr/bin/env python3
"""Check the periodic-window lower bounds in windows.json, using only Python's standard library.

    python3 verify_windows.py

Each entry gives a partition of Z_v into k sets that are Sidon modulo v -- for Singer,
the q+1 lines L - l (l in L) through the point 0, where L is the stored line (a perfect
difference set containing 0); for Ruzsa, the classes j = x - g^x (mod p) on
Z_{p-1} x Z_p = Z_{p(p-1)} -- together with an affine map x -> u*x + t and a length m.
The script checks that the classes are Sidon modulo v and partition Z_v, builds the
colouring of [v+m] given by x -> class of u*x + t (mod v), and checks every class of that
colouring by distinct pair sums and by distinct positive differences. It then prints
SR(k) >= v + m + 1.
"""
import itertools
import json
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent


def classes_mod(e):
    v, k = e["v"], e["k"]
    col = [None] * v
    if e["construction"] == "Singer":
        L = e["line"]
        assert len(L) == k and 0 in L
        diffs = sorted((a - b) % v for a in L for b in L if a != b)
        assert diffs == list(range(1, v)), "line is not a perfect difference set"
        for j, l in enumerate(L):
            for a in L:
                r = (a - l) % v
                if r:
                    assert col[r] is None
                    col[r] = j
        col[0] = -1                     # the common point: on every line through 0
    else:
        p, g = e["p"], e["g"]
        assert v == p * (p - 1) and k == p
        assert len({pow(g, i, p) for i in range(p - 1)}) == p - 1, "g is not a primitive root"
        col = [((x % p) - pow(g, x % (p - 1), p)) % p for x in range(v)]
    for c in range(k):                  # Sidon modulo v (the common point may join any class)
        cls = [x for x in range(v) if col[x] in (c, -1)]
        s = [(a + b) % v for a, b in itertools.combinations_with_replacement(cls, 2)]
        assert len(s) == len(set(s)), "class not Sidon modulo v"
    return col


def main():
    for e in json.loads((HERE / "windows.json").read_text()):
        v, k, m, u, t = e["v"], e["k"], e["m"], e["u"], e["t"]
        assert gcd(u, v) == 1 and 0 <= m <= k
        col = classes_mod(e)
        N = v + m
        res = [col[(u * x + t) % v] for x in range(1, N + 1)]
        used = {c for c in res[:m] if c >= 0}
        spare = next((c for c in range(k) if c not in used), 0)
        res = [spare if c < 0 else c for c in res]
        for c in range(k):
            A = [x for x in range(1, N + 1) if res[x - 1] == c]
            sums = [a + b for a, b in itertools.combinations_with_replacement(A, 2)]
            difs = [b - a for a, b in itertools.combinations(A, 2)]
            assert len(sums) == len(set(sums)) and len(difs) == len(set(difs)), "not a Sidon partition"
        assert N + 1 == e["lower_bound"]
        print(f"VERIFIED SR({k}) >= {N + 1}  ({e['construction']}, v={v}, m={m}, u={u}, t={t})", flush=True)


if __name__ == "__main__":
    main()
