# Verification code for "Bounds for Sidon-Ramsey and van der Waerden numbers"

Edgardo Roldán, Miguel Raggi and Tania Garibay.

This repository contains what a referee needs to check the explicit examples and finite computations of the paper, and, in `search/`, the programs of the exact search of Section 5.4, so that its runs (including the negative results) can be repeated. The checks use Python's standard library (3.9 or later); the upper-bound checks also need a C++17 compiler with OpenMP. All arithmetic in the checks is exact (integers or rationals), and no solver, SAT or floating-point computation is involved.

The verifications have not been formalized in a proof assistant.

## What is checked, and how to run it

| Command | What it checks | Time |
|---|---|---|
| `python3 verify.py --lower-only` | Each partition in `lower_partitions.json` (the Sidon-Ramsey lower bounds SR(8) to SR(20) of Table 1): it partitions `[1, n]` into at most `k` parts, each part has all pair sums distinct and all positive differences distinct, and its SHA-256 matches the stored one. | seconds |
| `python3 verify_windows.py` | The periodic-window bounds in `windows.json` (Singer and Ruzsa constructions): the classes are Sidon modulo `v` and partition `Z_v`, and the colouring of `[v+m]` they induce has Sidon classes. Prints SR(k) >= v+m+1 for every row. | under a second |
| `python3 check_windows_max.py` | For every row of `windows.json`, that `m` is the longest window over all affine images `x -> u*x + t`, and that it is at least the length `m_0` guaranteed by the averaging bound, Proposition 4.3 (the columns of the table of windows that `verify_windows.py` does not check). | about 40 seconds |
| `python3 vdw/verify.py` | The van der Waerden colourings in `vdw/certificates.json` have no monochromatic 3-term progression; the two covers of `Z_19` and `Z_37`; their products with Heule's colourings of `[170]` and `[75]` are the stored colourings of `[3230]` and `[2775]`; the facts about `Z_19` and `Z_27` used in the paper; and `r_3(Z_q)` for `q <= 20`. Ends with `ALL OK`. | under a second |
| `python3 verify.py --threads N` | Everything in the first row, plus the upper bound SR(9) <= 109: regenerates the complete lists of 13- and 14-element Sidon subsets of `[109]` and the four lists of 12-element subsets needed in the proof (C++ program `enumerate_sums.cpp`), and checks the integer weights in `certificates.json` against every set of the family, the disjoint pairs and the reflection identities. Intermediate files go in `verification/` (not tracked). | about two minutes on a 24-core machine |
| `python3 check_tiling108.py` | The fractional tiling of `[108]` by 108 twelve-element Sidon sets with rational weights (`tiling108.json`), which shows why the weight method stops at 109 (end of Section 3.3). | under a second |
| `python3 verify_upper.py --threads N` | The upper bounds SR(10) <= 132, SR(11) <= 158 and SR(12) <= 185 (Theorem 1.2, Section 3.4). Lists the Golomb rulers the proof needs with the C++ program `rulers.cpp`, which first proves the optimal lengths g_3, ..., g_12 it uses to prune and gets g_13, ..., g_16 as the shortest rulers of its lists. Checks that no four of the 15-element Sidon subsets of `[158]` and no five of the 16-element ones of `[185]` are pairwise disjoint. Checks the 142 signed-weight certificates in `certificates_sr10.json` against every 13-element Sidon set inside each complement (C++ program `sets_inside.cpp`). With `--g17` it also shows that no 17-mark ruler has length at most 184 (otherwise g_17 = 199 is taken from the literature); with `--exhaustive` it also searches the 142 complements exhaustively (`exact_cover.cpp`), after a positive control on the partition of `[102]`. | about 10 minutes on a 24-core Intel Core i9-14900KF; 11 with `--g17`, about 30 with `--exhaustive` |
| `g++ -O2 -std=c++17 vdw/r3_cyclic.cpp -o r3 && ./r3 lo hi` | Exhaustive computation of `r_3(Z_q)` for `lo <= q <= hi` (the values for `21 <= q <= 44` are stored in `vdw/r3_cyclic.json`). | seconds |

## Files

- `verify.py`, `enumerate_sums.cpp`: the upper bound SR(9) <= 109 and the lower-bound partitions.
- `certificates.json`: integer weights (with their common denominator) and the six 13-element sets (`large_sets`) used in the upper-bound proof.
- `lower_partitions.json`: the strongest verified partition for each number of colours.
- `windows.json`, `verify_windows.py`, `check_windows_max.py`: periodic-window bounds for more colours, and the maximality of each window.
- `tiling108.json`, `check_tiling108.py`: the fractional tiling of `[108]`.
- `verify_upper.py`, `rulers.cpp`, `sets_inside.cpp`, `exact_cover.cpp`, `certificates_sr10.json`: the upper bounds for ten, eleven and twelve colours.
- `search/`: the exact search for `[81]` and `[102]`, with build commands, positive controls and the commands of the runs (see `search/README.md`).
- `vdw/`: van der Waerden colourings (`certificates.json`), the covers of `Z_19` and `Z_37`, Heule's certificates (`heule_W_5_3_170.cert`, `heule_W_4_3_37.cert`) and the `r_3(Z_q)` computation.

## Notes

- The partitions are witnesses: they can be checked regardless of whether the randomized searches that found them can be replayed. The fields `checkpoint` and `parent` in `lower_partitions.json` give the names of the files they came from in the authors' search runs (directories removed); they are provenance only and the verifier does not read them.
- Heule's colourings come from his repository, cited in the paper.

## Reference machine

The two-minute figure for `verify.py` is for a 24-core Intel Core Ultra 9 285K; the figures for `verify_upper.py` are for a 24-core Intel Core i9-14900KF.

## Acknowledgments

This work was supported by UNAM-PAPIIT project IN114726.

## License

MIT (see `LICENSE`), covering the code and the data files. The colourings of Heule (`vdw/heule_*.cert`) are his, from the repository cited in the paper.
