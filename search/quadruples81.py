#!/usr/bin/env python3
"""Print one start line for completar.c per family of four disjoint 11-element Sidon sets in [81].

    python3 quadruples81.py ss11_81.txt > q81.txt     # 12,918 lines

A start line has 81 comma-separated entries: -c puts the number in the fixed class c
(c = 1..4, the four 11-element sets), 0 leaves it free for the four open classes.
"""
import sys

ss = [frozenset(map(int, line.split(","))) for line in open(sys.argv[1])]
bits = [sum(1 << x for x in s) for s in ss]
later = [[j for j in range(i + 1, len(ss)) if not bits[i] & bits[j]] for i in range(len(ss))]
count = 0
for i in range(len(ss)):
    for j in later[i]:
        u = bits[i] | bits[j]
        for k in later[j]:
            if u & bits[k]:
                continue
            v = u | bits[k]
            for l in later[k]:
                if v & bits[l]:
                    continue
                col = [0] * 81
                for c, t in enumerate((i, j, k, l)):
                    for x in ss[t]:
                        col[x - 1] = -(c + 1)
                print(",".join(map(str, col)))
                count += 1
print(count, "quadruples", file=sys.stderr)
