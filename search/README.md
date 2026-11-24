# The exact search (Section 5.4 of the paper)

These are the programs of the exact search that found the partitions of `[81]` into eight
Sidon sets and of `[102]` into nine, and that gave the negative results of Section 5.4. The
partitions themselves are certified by `../verify.py --lower-only` and do not depend on this
code; it is here so that the searches, including the negative results, can be repeated.

| File | What it does |
|---|---|
| `ssg.c` | Lists all K-element Sidon subsets of `[N]` into `ssK_N.bin` (K bytes per set). Prunes with the optimal Golomb ruler lengths; OpenMP. |
| `tuplas3.c` | Lists the R-tuples of pairwise disjoint sets from a list (one set per line, comma separated), up to the reflection `x -> N+1-x`: a tuple is kept if its sorted index tuple is at most the sorted tuple of the reflected indices. It can stream to stdout in `P` shards. |
| `mkl_stream.py` | Turns a line of indices `i_j_k` into the line `i_j_k a1,a2,...` that `xc9` reads. |
| `xc9.c` | For each line (a tuple of fixed largest classes), covers the rest of `[N]` with sets of the next size (Knuth's Algorithm X, branching on the point with the fewest options), while the points of the smaller classes go to one leftover pool that must split into Sidon sets of prescribed sizes. Prints `FOUND`, `EXHAUSTED` or `TIMEOUT` per profile, and the first solution as a `SOL` line on stderr. Candidate size `-DSZ` (default 11) and fixed-class size `-DGS` (default 12). |
| `completar.c`, `quadruples81.py` | Exact completion with some classes fixed, colouring one free point at a time; used for the four 11-classes of `[81]`. The script writes its start lines. |
| `disjoint.c` | Counts the families of R pairwise disjoint sets in a list. It shares no code with the programs above; it decides the points in increasing order (uncovered, or covered by a set starting there). |

The profile argument of `xc9` is a list of profiles separated by `/`; `4:8/3:9,9` means: four
classes of size `SZ` plus a leftover class of size 8, or three of size `SZ` plus two of size 9.

## Build

    gcc -O3 -march=native -fopenmp -o ssg ssg.c
    gcc -O3 -march=native -o tuplas3 tuplas3.c
    gcc -O3 -march=native -o xc9 xc9.c -lm                    # [102]: SZ=11, GS=12
    gcc -O3 -march=native -DSZ=10 -DGS=11 -o xc8 xc9.c -lm    # [81]:  SZ=10, GS=11
    gcc -O3 -march=native -fopenmp -o disjoint disjoint.c
    # one set per line, for tuplas3 and mkl_stream.py:
    tolines() { python3 -c "import sys; k=int(sys.argv[2]); b=open(sys.argv[1],'rb').read(); print(''.join(','.join(map(str,b[i:i+k]))+'\n' for i in range(0,len(b),k)),end='')" "$1" "$2"; }

## Set lists

    ./ssg 81 11; ./ssg 81 10        # 4,030 and 11,817,176 sets
    ./ssg 102 12; ./ssg 102 11      # 43,698 and 125,628,796 sets (ss11_102.bin is 1.4 GB)
    tolines ss11_81.bin 11 > ss11_81.txt
    tolines ss12_102.bin 12 > ss12_102.txt

## Positive controls (seconds)

The partitions in `../lower_partitions.json` are found again from their largest classes:

    # [81]: fix its three 11-classes; FOUND with profile 3:9,9 (0.6 s, exactly one solution)
    echo "c81 <the 33 elements of the three 11-classes, comma separated>" | ./xc8 81 ss10_81.bin batch "4:8/3:9,9" 60
    # [102]: fix its four 12-classes; FOUND with profile 4:10 (about 7 s)
    echo "c102 <the 48 elements of the four 12-classes, comma separated>" | ./xc9 102 ss11_102.bin batch "4:10" 120

and `./disjoint ss11_81.bin 81 11 4` prints 12918, the number of families of four disjoint
11-element Sidon sets in `[81]` (6,483 up to reflection, as `./tuplas3 81 ss11_81.txt 4 out` reports).

## The runs of the paper

**[81], eight colours, three 11-classes (exhaustive).** 621,365 canonical triples
(`./tuplas3 81 ss11_81.txt 3 t3.txt`), both profiles of the five remaining classes:

    python3 mkl_stream.py ss11_81.txt < t3.txt | ./xc8 81 ss10_81.bin batch "4:8/3:9,9" 600

The run of the paper split the triples into 16 processes and took 19.4 CPU-hours; exactly five
triples are `FOUND`, each with profile 3:9,9 and one solution.

**[81], four or more 11-classes (exhaustive).** This case was decided by a different program,
`completar.c`, which colours the free points one at a time (the one with the fewest classes
left first) and so needs no list of 10-element sets. All 12,918 families of four disjoint
11-element Sidon sets are run, not only one per reflection pair:

    gcc -O2 -march=native -o completar completar.c
    python3 quadruples81.py ss11_81.txt > q81.txt                  # 12,918 start lines
    split -n l/60 -d q81.txt part_ && for f in part_??; do ./completar 81 8 1 < $f > $f.out & done; wait
    grep -c 'completaciones=0 ' part_??.out                          # adds up to 12918

That took 272 CPU-hours (median 69 s per quadruple). Positive control: with the three
11-classes and one 10-class of the partition in `../lower_partitions.json` fixed (start line
with `-1`..`-4` on their points and `0` elsewhere), `./completar 81 8 1` prints a `SOL` line and
`completaciones=1+` (about 95 s).

**[102], nine colours.** With t fixed 12-classes, `12t + 11(9-t) >= 102` gives `t >= 3`. The
tuples are streamed in 14 shards (`r = 0..13`):

    ./tuplas3 102 ss12_102.txt 4 - 14 r | python3 mkl_stream.py ss12_102.txt \
        | ./xc9 102 ss11_102.bin batch "4:10" 60 > out4_r 2> err4_r          # 4-tuples
    ./tuplas3 102 ss12_102.txt 5 - 14 r | python3 mkl_stream.py ss12_102.txt \
        | ./xc9 102 ss11_102.bin batch "3:9/2:10,10" 60                      # 5-tuples
    ./tuplas3 102 ss12_102.txt 6 - 14 r | python3 mkl_stream.py ss12_102.txt \
        | ./xc9 102 ss11_102.bin batch "2:8/1:10,9/0:10,10,10" 60            # 6-tuples

All 11,366,148 five-tuples and all 286 six-tuples are `EXHAUSTED`. The four-tuple run (about
1.5e9 cases, estimated) found the partition of the paper after about 0.9e9 cases and six days
on 14 processes. Seven or more 12-classes are impossible:

    ./disjoint ss12_102.bin 102 12 7      # prints 0 (about 3 minutes on 30 threads)
