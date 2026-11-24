// Exact completion of a partition whose large classes are fixed: can the free numbers
// of [1,N] be split among the open classes so that every class is a Sidon set?
//
// This is the second half of the method of Espinosa-García and Pellicer
// (arXiv:2309.08553): enumerate tuples of disjoint large Sidon sets, then decide each
// one by exhaustive search. repair.c with fixed classes does the same job heuristically,
// so a failure there proves nothing; a failure here is a proof for that tuple.
//
// The search colours one free number at a time, always the one with the fewest classes
// left (a number with none ends the branch). A class stays Sidon while its pairwise
// differences are distinct, so each class keeps its difference set as a bitmask, and
// its members twice: as a mask and mirrored (bit N+1-a), which turns the differences
// of x with every member into two shifts.
// The open classes start empty and are interchangeable, so a number may only open the
// first empty one; every partition is then found exactly once.
//
// With a cap, no open class may grow past cap elements. That is how the cases "exactly
// a classes of the largest size" are searched: fix those a classes, cap the rest one
// below. The cap also prunes. Class c can still take at most min(cap - |c|, the free
// numbers that fit c), and if these do not add up to the numbers still free, the
// branch is dead. With no slack, as in [81] with one 11-class fixed and seven classes
// capped at 10 for 70 numbers, that bites at every node.
//
//   completar N K [limit [cap]] < starts
//
// Each line of starts is one start file in repair.c's format: N comma-separated
// entries, -c fixing that number in class c (closed to everyone else), 0 leaving it
// free. The fixed classes must be Sidon sets. For each line it prints the number of
// completions (stopping at limit, default 1; 0 counts them all) and the search nodes.
// Each completion found is printed as a colouring line beginning with "SOL".
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#define MAXN 120
#define MAXK 16
typedef unsigned __int128 u128;

static int N, K, nopen, nfree;
static int col[MAXN + 1], freex[MAXN];
static int mem[MAXK][MAXN], msz[MAXK];
static u128 diffs[MAXK], memb[MAXK], mirr[MAXK];
static long long nodes, found, limit;
static int cap = MAXN;

// The differences x would add to class c, or 0 with *ok = 0 if they clash. x - a and
// b - x coincide when a + b = 2x, so the new differences are checked among themselves too.
static inline u128 adds(int x, int c, int *ok) {
    u128 lo = mirr[c] >> (N + 1 - x), hi = memb[c] >> x;   // {x-a : a<x}, {a-x : a>x}
    *ok = !(lo & hi) && !((lo | hi) & diffs[c]);
    return lo | hi;
}

static inline void put(int x, int c, u128 d) {
    diffs[c] |= d; memb[c] |= (u128)1 << x; mirr[c] |= (u128)1 << (N + 1 - x);
    mem[c][msz[c]++] = x; col[x] = c;
}

static inline void take(int x, int c, u128 d) {
    diffs[c] &= ~d; memb[c] &= ~((u128)1 << x); mirr[c] &= ~((u128)1 << (N + 1 - x));
    msz[c]--; col[x] = -1;
}

static void report(void) {
    printf("SOL ");
    for (int x = 1; x <= N; x++) printf("%d%s", col[x] + 1, x < N ? "," : "\n");
}

static void search(void) {
    if (limit && found >= limit) return;
    nodes++;
    // most constrained free number; classes are tried in index order
    int best = -1, bestn = MAXK + 1, opts[MAXK], nb = 0;
    int first_empty = -1, left = 0, fits[MAXK] = {0};
    for (int c = K - nopen; c < K; c++) if (!msz[c]) { first_empty = c; break; }
    for (int i = 0; i < nfree; i++) {
        int x = freex[i];
        if (col[x] >= 0) continue;
        left++;
        int o[MAXK], n = 0, ok;
        for (int c = K - nopen; c < K; c++) {
            if (msz[c] >= cap || (!msz[c] && c != first_empty)) continue;
            adds(x, c, &ok);
            if (ok) { o[n++] = c; fits[c]++; }
        }
        if (n < bestn) { bestn = n; best = x; nb = n; memcpy(opts, o, n * sizeof *o); if (!n) return; }
    }
    if (best < 0) { found++; report(); return; }
    if (cap < MAXN) {
        int room = 0;
        for (int c = K - nopen; c < K; c++) {
            int r = cap - msz[c], f = msz[c] ? fits[c] : left;   // every number fits an empty class
            room += r < f ? r : f;
        }
        if (room < left) return;
    }
    for (int j = 0; j < nb; j++) {
        int c = opts[j], ok;
        u128 d = adds(best, c, &ok);
        put(best, c, d);
        search();
        take(best, c, d);
        if (limit && found >= limit) return;
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "uso: %s N K [limite [tope]] < inicios\n", argv[0]); return 2; }
    N = atoi(argv[1]); K = atoi(argv[2]); limit = argc > 3 ? atoll(argv[3]) : 1;
    if (argc > 4 && atoi(argv[4]) > 0) cap = atoi(argv[4]);
    char line[8192];
    while (fgets(line, sizeof line, stdin)) {
        int given[MAXN + 1], n = 0, nfix = 0;
        char *p = line;
        while (n < N) {
            char *e; long v = strtol(p, &e, 10);
            if (e == p) break;
            given[++n] = (int)v; p = e; while (*p == ',' || *p == ' ') p++;
        }
        if (n != N) { fprintf(stderr, "linea con %d entradas, se esperaban %d\n", n, N); return 2; }
        memset(msz, 0, sizeof msz); memset(diffs, 0, sizeof diffs);
        memset(memb, 0, sizeof memb); memset(mirr, 0, sizeof mirr);
        // fixed classes are renumbered to 0.. in order of appearance
        int map[MAXK + 1]; for (int c = 0; c <= MAXK; c++) map[c] = -1;
        nfree = 0;
        for (int x = 1; x <= N; x++) {
            col[x] = -1;
            if (given[x] > 0) { fprintf(stderr, "solo se aceptan -c y 0\n"); return 2; }
            if (!given[x]) { freex[nfree++] = x; continue; }
            int c = -given[x];
            if (map[c] < 0) map[c] = nfix++;
            int ok; u128 d = adds(x, map[c], &ok);
            if (!ok) { fprintf(stderr, "la clase fija %d no es de Sidon\n", c); return 2; }
            put(x, map[c], d);
        }
        nopen = K - nfix;
        nodes = found = 0;
        clock_t t0 = clock();
        search();
        printf("completaciones=%lld%s nodos=%lld segs=%.3f\n", found, limit && found >= limit ? "+" : "",
               nodes, (double)(clock() - t0) / CLOCKS_PER_SEC);
        fflush(stdout);
    }
    return 0;
}
