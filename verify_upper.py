#!/usr/bin/env python3
"""Check the upper bounds SR(10) <= 132, SR(11) <= 158 and SR(12) <= 185 (Section 3.4 of the paper).

    python3 verify_upper.py [--threads N] [--g17]

Needs a C++17 compiler with OpenMP and Python's standard library. Intermediate files go in
verification/.

1. Golomb rulers. rulers.cpp lists every M-mark Golomb ruler of length at most L. It cuts
   branches only with lower bounds for rulers with fewer marks, and every bound it uses was
   proved by an earlier run of this script: for M = 3, ..., 12 a run with L = g_M shows that
   the shortest M-mark ruler has the published length g_M. Then the script lists the 13- and
   14-mark rulers of length at most 131, the 15-mark rulers of length at most 157 and the
   16-mark rulers of length at most 184; their shortest members give g_13, ..., g_16. With
   --g17 it also shows that no 17-mark ruler has length at most 184 (about a minute);
   without it, that fact is taken from g_17 = 199 in the literature.
2. SR(11) <= 158 and SR(12) <= 185. Every Sidon subset of [158] has at most 15 elements, and
   every one of [185] at most 16, so a partition into 11 (12) Sidon sets has at least 4 (5)
   pairwise disjoint classes of that size. The script lists these sets and shows that no
   4 (5) of them are pairwise disjoint.
3. SR(10) <= 132. A partition into 10 Sidon sets has at least 2 classes with 14 elements.
   For each family of 2 or 3 pairwise disjoint 14-element Sidon subsets of [132] (there is
   none of 4), the rest of the partition would be 8 classes with 13 elements covering the
   complement U exactly (2 classes of 14), or 6 classes with 13 elements and one with 12
   (3 classes of 14). The signed weights w and the integer Z in certificates_sr10.json rule
   this out: every 13-element Sidon set inside U has weight w(S) >= 0 (checked by
   sets_inside.cpp against the complete list), and w(U) is smaller than the least weight that
   f = |U| - 13 * (number of 13-element classes) further points of U can have, namely
   f*Z - sum_{x in U} max(0, Z - w(x)); for f = 0 the condition is w(U) < 0.
4. With --exhaustive, a second check of SR(10) <= 132 that does not use the weights:
   exact_cover.cpp searches every one of the 142 complements exhaustively and must find no
   decomposition. Positive control first: with the four 12-element classes of the partition of
   [102] in lower_partitions.json fixed, it must find exactly that partition. The largest
   complements take about 20 minutes each on one core; they run in parallel.
"""
import argparse
import concurrent.futures
import itertools
import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "verification"
PUBLISHED = {2: 1, 3: 3, 4: 6, 5: 11, 6: 17, 7: 25, 8: 34, 9: 44, 10: 55, 11: 72, 12: 85}


def log(msg):
    print(msg, flush=True)


def compile_program(name):
    exe = OUT / name
    subprocess.run(["g++", "-O3", "-march=native", "-std=c++17", "-fopenmp", "-o", str(exe), str(HERE / f"{name}.cpp")], check=True)
    return exe


def is_golomb(r):
    d = [b - a for a, b in itertools.combinations(r, 2)]
    return len(set(d)) == len(d)


def list_rulers(exe, M, L, bounds, env):
    """All M-mark Golomb rulers of length <= L (both orientations), as tuples."""
    out = OUT / f"rulers_m{M}_L{L}.bin"
    t0 = time.time()
    res = subprocess.run([str(exe), str(M), str(L), ",".join(str(bounds[r]) for r in range(1, M)), str(out)], check=True, capture_output=True, text=True, env=env)
    info = json.loads(res.stdout)
    data = out.read_bytes()
    rulers = [tuple(data[i:i + M]) for i in range(0, len(data), M)]
    # what the C++ program returned must be what it claims: sorted, distinct, Golomb, first mark 0, length <= L
    assert len(rulers) == info["count"] and rulers == sorted(set(rulers))
    assert all(r[0] == 0 and r[-1] <= L and is_golomb(r) for r in rulers)
    log(f"  {M}-mark rulers of length <= {L}: {len(rulers)} ({time.time() - t0:.0f} s)")
    return rulers


def translates(rulers, n):
    """The Sidon subsets of [n] that are translates of the given rulers."""
    return sorted({tuple(a + s for a in r) for r in rulers if r[-1] <= n - 1 for s in range(1, n - r[-1] + 1)})


def disjoint_families(sets, t):
    """The families of t pairwise disjoint sets, as sorted index tuples."""
    masks = [sum(1 << x for x in s) for s in sets]
    found = []

    def rec(start, used, chosen):
        if len(chosen) == t:
            found.append(tuple(chosen))
            return
        for i in range(start, len(masks)):
            if not masks[i] & used:
                rec(i + 1, used | masks[i], chosen + [i])

    rec(0, 0, [])
    return found


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--threads", type=int, default=os.cpu_count())
    ap.add_argument("--g17", action="store_true", help="also show that no 17-mark Golomb ruler has length <= 184")
    ap.add_argument("--exhaustive", action="store_true", help="also search the 142 complements for SR(10) exhaustively")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads))
    rulers_exe = compile_program("rulers")
    inside_exe = compile_program("sets_inside")

    log("1. Golomb rulers")
    g = {1: 0, 2: 1}  # one mark has length 0; two marks need length >= 1
    for M in range(3, 13):
        rl = list_rulers(rulers_exe, M, PUBLISHED[M], g, env)
        assert rl and min(r[-1] for r in rl) == PUBLISHED[M]
        g[M] = PUBLISHED[M]
    catalogue = {}
    for M, L in [(13, 131), (14, 131), (15, 157), (16, 184)]:
        catalogue[M] = list_rulers(rulers_exe, M, L, g, env)
        g[M] = min(r[-1] for r in catalogue[M])
    log("  shortest rulers: " + ", ".join(f"g_{M} = {g[M]}" for M in range(3, 17)))
    assert [g[M] for M in (13, 14, 15, 16)] == [106, 127, 151, 177]
    if args.g17:
        assert not list_rulers(rulers_exe, 17, 184, g, env)
        log("  no 17-mark ruler of length <= 184, so g_17 > 184")
    else:
        log("  g_17 = 199 > 184 is taken from the literature (run with --g17 to check g_17 > 184)")

    log("2. SR(11) <= 158 and SR(12) <= 185")
    for k, n, M in [(11, 158, 15), (12, 185, 16)]:
        need = n - k * (M - 1)  # a k-partition of [n] into sets of at most M elements has >= need classes of M
        sets = translates(catalogue[M], n)
        fams = disjoint_families(sets, need - 1)
        assert not disjoint_families(sets, need)
        log(f"  [{n}]: every Sidon subset has <= {M} elements; {len(sets)} with {M} elements; {len(fams)} families of {need - 1} pairwise disjoint ones and none of {need}, which a {k}-partition would need: SR({k}) <= {n}")

    log("3. SR(10) <= 132")
    n, k, M = 132, 10, 14
    sets = translates(catalogue[14], n)
    fam = {t: disjoint_families(sets, t) for t in (2, 3, 4)}
    assert not fam[4]
    log(f"  [{n}]: {len(sets)} Sidon subsets with 14 elements; {len(fam[2])} disjoint pairs, {len(fam[3])} disjoint triples, no 4 pairwise disjoint")
    cert = json.load(open(HERE / "certificates_sr10.json"))
    assert (cert["k"], cert["n"], cert["M"]) == (k, n, M)
    bykey = {(c["t"], tuple(sorted(tuple(s) for s in c["family"]))): c for c in cert["certificates"]}
    cases = [(t, c) for t in (2, 3) for c in fam[t]]
    assert len(bykey) == len(cert["certificates"]) == len(cases)
    lines, checks = [], []
    for t, c in cases:
        cw = bykey[(t, tuple(sorted(sets[i] for i in c)))]
        U = sorted(set(range(1, n + 1)) - set().union(*(sets[i] for i in c)))
        w = {int(x): int(v) for x, v in cw["weights"].items()}
        assert sorted(w) == U
        lines.append(" ".join(str(w[x]) if x in w else "x" for x in range(1, n + 1)))
        # the other 10 - t classes have <= 13 elements and cover U, which is 13*(10-t) - deficit points:
        # with deficit 0 all have 13 elements, with deficit 1 one class has 12 and the others 13
        deficit = 13 * (k - t) - len(U)
        assert deficit == (0 if t == 2 else 1) and len(U) == (104 if t == 2 else 90)
        free = len(U) - 13 * (len(U) // 13)  # points not in a 13-element class: 0 or 12
        checks.append((t, U, w, free, int(cw.get("Z", 0)), cw["sets_inside"]))
    (OUT / "sr10_cases.txt").write_text("\n".join(lines) + "\n")
    t0 = time.time()
    res = subprocess.run([str(inside_exe), str(OUT / "rulers_m13_L131.bin"), "13", str(n), str(OUT / "sr10_cases.txt")], check=True, capture_output=True, text=True, env=env)
    rows = [tuple(map(int, l.split())) for l in res.stdout.split("\n") if l.strip()]
    assert len(rows) == len(checks)
    sizes = {2: [], 3: []}
    for (t, U, w, free, Z, expected), (count, minw) in zip(checks, rows):
        assert count == expected and count > 0
        total = sum(w[x] for x in U)
        bound = 0 if free == 0 else free * Z - sum(max(0, Z - w[x]) for x in U)
        assert minw >= 0 and total < bound, (t, U)
        sizes[t].append(count)
    log(f"  {len(checks)} complements ({time.time() - t0:.0f} s): 104 points with {min(sizes[2])}-{max(sizes[2])} thirteen-element Sidon subsets (2 classes of 14, exact cover by 8), 90 points with {min(sizes[3])}-{max(sizes[3])} (3 classes of 14, 6 sets and 12 more points); every certificate holds: SR(10) <= 132")
    if args.exhaustive:
        exhaustive(args, env, g, inside_exe, [(t, U) for (t, U, *_rest) in checks])
    log("All checks passed.")


def exact_cover(exe, U, sets_file, M, free):
    """Number of ways to split U into members of sets_file plus `free` points forming a Sidon set."""
    ufile = Path(str(sets_file)[:-4] + "_U.txt")
    ufile.write_text(" ".join(map(str, U)) + "\n")
    out = subprocess.run([str(exe), str(ufile), str(sets_file), str(M), str(free)], check=True, capture_output=True, text=True).stdout
    head = out.split("\n")[0].split()
    solution = sorted(tuple(sorted(map(int, l.split()[1:]))) for l in out.split("\n") if l.startswith(("SET", "LEFT")))
    return int(head[head.index("solutions") + 1]), solution


def exhaustive(args, env, g, inside_exe, cases):
    log("4. SR(10) <= 132 again, by exhaustive search")
    cover_exe = compile_program("exact_cover")
    # positive control on the partition of [102]: fix its four 12-element classes
    parts = [sorted(c) for c in json.load(open(HERE / "lower_partitions.json"))["9"]["parts"]]
    fixed = [c for c in parts if len(c) == 12]
    assert len(fixed) == 4 and sorted(map(len, parts)) == [10, 11, 11, 11, 11, 12, 12, 12, 12]
    U = sorted(set(range(1, 103)) - set().union(*map(set, fixed)))
    r11 = list_rulers(OUT / "rulers", 11, 101, g, env)
    ctl = OUT / "control102"
    ctl.mkdir(exist_ok=True)
    (ctl / "case.txt").write_text(" ".join("0" if x in U else "x" for x in range(1, 103)) + "\n")
    subprocess.run([str(inside_exe), str(OUT / "rulers_m11_L101.bin"), "11", "102", str(ctl / "case.txt"), str(ctl)], check=True, capture_output=True, env=env)
    count, solution = exact_cover(cover_exe, U, ctl / "0.bin", 11, 10)
    assert count == 1 and solution == sorted(tuple(c) for c in parts if len(c) < 12), (count, solution)
    log(f"  control: [102] with its four 12-element classes fixed has exactly one completion, the stored partition ({len(r11)} rulers with 11 marks)")
    sets_dir = OUT / "sr10_sets"
    sets_dir.mkdir(exist_ok=True)
    subprocess.run([str(inside_exe), str(OUT / "rulers_m13_L131.bin"), "13", "132", str(OUT / "sr10_cases.txt"), str(sets_dir)], check=True, capture_output=True, env=env)
    t0 = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as pool:
        jobs = [pool.submit(exact_cover, cover_exe, U, sets_dir / f"{i}.bin", 13, 0 if t == 2 else 12) for i, (t, U) in enumerate(cases)]
        counts = [j.result()[0] for j in jobs]
    assert counts == [0] * len(cases), counts
    log(f"  none of the {len(cases)} complements splits as required ({time.time() - t0:.0f} s): SR(10) <= 132")


if __name__ == "__main__":
    sys.exit(main())
