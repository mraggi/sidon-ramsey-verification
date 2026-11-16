"""Check every claim of the paper that is a finite computation. Python standard library only.

    python3 vdw/verify.py       # from the repository root

1. Each colouring in vdw/certificates.json uses at most r colours and has no
   monochromatic 3-term progression: for each colour class C (sorted), every pair x < z in C
   with x + z even and (x + z)/2 in C would be one; there must be none.
2. The covers in vdw/cover_z19.json and cover_z37.json: sixteen subsets of Z_19,
   resp. fifteen of Z_37, none containing a 3-term progression x, x+d, x+2d (d != 0 mod q),
   together covering every residue at least five, resp. four, times.
3. The product of each cover with a base colouring (Kim's palette theorem) is exactly the
   stored product colouring: Heule's W(5,3) > 170 certificate gives W(16,3) > 3230, and
   Heule's 74-point 4-colouring (file W_4_3_37.cert) with one point of colour b put in front,
   a valid 4-colouring of 75 points, gives W(15,3) > 2775.
4. The facts about Z_19 stated in the paper: it has exactly 57 progression-free 6-subsets,
   none of size 7, and the 57 form one orbit of {0,1,3,12,14,15} under x -> a x + b.
5. Every progression-free 8-subset of Z_27 has four elements in each of two classes mod 3 and
   none in the third (used to show that no 17 such sets cover Z_27 five times).
6. r_3(Z_q) for q <= 20, by exhaustive search (the larger values in r3_cyclic.json come from
   the C++ program r3_cyclic.cpp, which is the same search with bit masks).
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def decode(word: str) -> list[int]:
    return [ord(ch) - ord("a") for ch in word]


def mono_3aps(col: list[int]) -> int:
    classes: dict[int, list[int]] = {}
    for x, c in enumerate(col):
        classes.setdefault(c, []).append(x)
    bad = 0
    for pos in classes.values():
        members = set(pos)
        for i, x in enumerate(pos):
            for z in pos[i + 1:]:
                if (x + z) % 2 == 0 and (x + z) // 2 in members:
                    bad += 1
    return bad


def cyclic_progressions(q: int):
    return {(x, (x + d) % q, (x + 2 * d) % q) for x in range(q) for d in range(1, q)}


def progression_free(A, q: int, progs=None) -> bool:
    S = set(A)
    return not any(a in S and b in S and c in S for a, b, c in (progs or cyclic_progressions(q)))


def main() -> int:
    ok = True
    data = json.load(open(HERE / "certificates.json"))
    for r, rec in data["best"].items():
        col = decode(rec["coloring"])
        good = len(col) == rec["n"] and max(col) < int(r) and mono_3aps(col) == 0
        print(f"W({r},3) > {rec['n']}: {'ok' if good else 'FAILED'}")
        ok &= good
    bases = {16: ("cover_z19.json", decode(open(HERE / "heule_W_5_3_170.cert").read().strip()), 170),
             15: ("cover_z37.json", [1] + decode(open(HERE / "heule_W_4_3_37.cert").read().strip()), 75)}
    for rr, (cover_file, beta, n_beta) in bases.items():
        prod = decode(data[f"product{rr}"]["coloring"])
        good = len(prod) == data[f"product{rr}"]["n"] and max(prod) < rr and mono_3aps(prod) == 0
        print(f"W({rr},3) > {len(prod)} (product): {'ok' if good else 'FAILED'}")
        ok &= good

        cover = json.load(open(HERE / cover_file))
        q, r, s, sets = cover["q"], cover["r"], cover["s"], cover["sets"]
        progs = cyclic_progressions(q)
        good = r == rr and len(sets) == r and all(progression_free(A, q, progs) for A in sets) \
            and all(sum(u in A for A in sets) >= s for u in range(q))
        print(f"cover of Z_{q} by {r} progression-free sets, every residue at least {s} times: {'ok' if good else 'FAILED'}")
        ok &= good

        good_beta = len(beta) == n_beta and max(beta) < s and mono_3aps(beta) == 0
        palettes = [[c for c in range(r) if u in sets[c]][:s] for u in range(q)]
        chi = [palettes[u][beta[t]] for t in range(len(beta)) for u in range(q)]
        good = good_beta and chi == prod
        print(f"base W({s},3) > {n_beta} colouring valid, and product = stored W({r},3) > {q * n_beta} colouring: "
              f"{'ok' if good else 'FAILED'}")
        ok &= good
    q = 19
    progs = cyclic_progressions(q)

    six = [A for A in itertools.combinations(range(q), 6) if progression_free(A, q, progs)]
    seven = sum(1 for A in itertools.combinations(range(q), 7) if progression_free(A, q, progs))
    base = (0, 1, 3, 12, 14, 15)
    orbit = {tuple(sorted((a * x + b) % q for x in base)) for a in range(1, q) for b in range(q)}
    good = len(six) == 57 and seven == 0 and set(six) == orbit
    print(f"Z_19: {len(six)} progression-free 6-sets, {seven} 7-sets, one affine orbit of {base}: {'ok' if good else 'FAILED'}")
    ok &= good

    q, inv2 = 27, 14                                   # 2 * 14 = 1 mod 27
    profiles, sizes = set(), [0] * 10

    def extend(cur, banned, start):
        sizes[len(cur)] += 1
        if len(cur) == 8:
            profiles.add(tuple(sorted((sum(1 for x in cur if x % 3 == i) for i in range(3)), reverse=True)))
        if len(cur) == 9:
            return
        for x in range(start, q):
            if banned[x] or len(cur) + (q - x) < 8:
                continue
            new = [(2 * y - x) % q for y in cur] + [(2 * x - y) % q for y in cur] + [(x + y) * inv2 % q for y in cur]
            for z in new:
                banned[z] += 1
            extend(cur + [x], banned, x + 1)
            for z in new:
                banned[z] -= 1
    extend([], [0] * q, 0)
    good = profiles == {(4, 4, 0)} and sizes[8] == 2187 and sizes[9] == 0
    print(f"Z_27: {sizes[8]} progression-free 8-sets, none of size 9, all 8-sets split 4+4+0 mod 3: {'ok' if good else 'FAILED'}")
    ok &= good

    r3 = {int(k): v for k, v in json.load(open(HERE / "r3_cyclic.json")).items()}
    for qq in range(2, 21):
        pr = cyclic_progressions(qq)
        size = r3[qq]
        has = any(progression_free(A, qq, pr) for A in itertools.combinations(range(qq), size))
        none_bigger = size == qq or not any(progression_free(A, qq, pr) for A in itertools.combinations(range(qq), size + 1))
        if not (has and none_bigger):
            print(f"r_3(Z_{qq}) = {size}: FAILED")
            ok = False
    print("r_3(Z_q) for q <= 20:", "ok" if ok else "see above")
    print("ALL OK" if ok else "SOME CHECK FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
