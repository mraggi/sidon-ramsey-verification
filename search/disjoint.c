// disjoint SETS.bin N K R : count the families of R pairwise disjoint sets in a list of
// K-element subsets of [N] (K bytes per set, as written by ssg).
//
// Independent of tuplas3 and xc9. It decides the points 1, 2, ..., N in order: the smallest
// undecided point x is either left uncovered, or covered by a set whose smallest element
// is x (a set that contains x and a smaller point y was already decided at y). At most
// N - R*K points can stay uncovered, which prunes the search. Every family is counted
// exactly once, as a set of sets (not up to reflection).
//
//   gcc -O3 -march=native -fopenmp -o disjoint disjoint.c
//   ./ssg 102 12 && ./disjoint ss12_102.bin 102 12 7    # 0
#include <stdio.h>
#include <stdlib.h>
typedef unsigned __int128 u128;
static int N, K, budget;
static u128 *S; static long nS; static long first[130];  // first[x] = first set whose minimum is >= x

static int minel(u128 m) { int x = 0; while (!((m >> x) & 1)) x++; return x; }
static int cmp(const void *a, const void *b) {
    int x = minel(*(const u128 *)a), y = minel(*(const u128 *)b); return (x > y) - (x < y);
}

static long rec(int x, u128 used, int skipped, int left) {
    if (left == 0) return 1;
    while (x <= N && ((used >> x) & 1)) x++;
    if (x > N) return 0;
    long c = 0;
    if (skipped < budget) c += rec(x + 1, used, skipped + 1, left);
    for (long i = first[x]; i < first[x + 1]; i++)
        if (!(S[i] & used)) c += rec(x + 1, used | S[i], skipped, left - 1);
    return c;
}

int main(int argc, char **argv) {
    if (argc != 5) { fprintf(stderr, "usage: %s SETS.bin N K R\n", argv[0]); return 2; }
    N = atoi(argv[2]); K = atoi(argv[3]); int R = atoi(argv[4]); budget = N - R * K;
    if (N > 127 || K > 64 || R < 1 || budget < 0) { fprintf(stderr, "need N <= 127, R >= 1, R*K <= N\n"); return 2; }
    FILE *f = fopen(argv[1], "rb"); if (!f) { perror(argv[1]); return 2; }
    long cap = 1 << 16; S = malloc(cap * sizeof *S); unsigned char b[64];
    while (fread(b, 1, K, f) == (size_t)K) {
        if (nS == cap) { cap *= 2; S = realloc(S, cap * sizeof *S); }
        u128 m = 0; for (int i = 0; i < K; i++) { if (b[i] < 1 || b[i] > N) { fprintf(stderr, "bad element\n"); return 2; } m |= (u128)1 << b[i]; }
        S[nS++] = m;
    }
    fclose(f);
    qsort(S, nS, sizeof *S, cmp);
    for (int x = 1, i = 0; x <= N + 1; x++) { while (i < nS && minel(S[i]) < x) i++; first[x] = i; }
    fprintf(stderr, "%ld sets of size %d in [%d]\n", nS, K, N);
    long total = 0;
    // one task per possible first set; the points before its minimum stay uncovered
    #pragma omp parallel for schedule(dynamic, 1) reduction(+ : total)
    for (long i = 0; i < nS; i++) {
        int x0 = minel(S[i]);
        if (x0 - 1 <= budget) total += rec(x0 + 1, S[i], x0 - 1, R - 1);
    }
    printf("%ld\n", total);
    return 0;
}
