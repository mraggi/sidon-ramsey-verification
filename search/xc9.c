// Algorithm X with MRV where only the 10-element classes are candidate sets; the small
// classes are a single "leftover" pool that must split into Sidon sets of given sizes.
//   xcover3 N SETS10.bin K10 SMALL MODE SECS < fixed
//     K10   = number of 10-classes to place, SMALL = sizes of the small classes ("8" or "9,9",
//     "-" for none); MODE 0 = first, 1 = count all, >=2 = Knuth estimator with MODE probes.
// At a node each uncovered number x has options: every remaining 10-set containing x, plus
// "x goes to the leftover" if the leftover with x can still split into the small classes.
// Branch on the number with fewest options. Numbers no 10-set can cover are forced to the
// leftover at once. When all K10 sets are placed the rest is the leftover; check it exactly.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
typedef unsigned __int128 u128;
#ifndef SZ
#define SZ 11
#endif
#ifndef GS
#define GS 12
#endif
static int N, MODE, K10, nsmall, small[8], smalltot;
static double deadline, t0;
static u128 *M; static long nM;
static int **lev; static long long nodes, sols, nodes_at[32], splitcalls;
static int timeout; static u128 chosen[16]; static int nch; static int printed;
static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+t.tv_nsec*1e-9; }
static inline int pc(u128 m){ return __builtin_popcountll((unsigned long long)m) + __builtin_popcountll((unsigned long long)(m >> 64)); }

// can the numbers of mask m be split into nsmall Sidon sets of sizes <= small[i]
// (exact sizes when exact=1)?  Plain backtracking; m has at most ~27 elements.
static int el[128], ne, col[40], csz[8]; static u128 cdif[8]; static int cmem[8][16];
static int sp(int i, int exact) {
    if (i == ne) { if (!exact) return 1; for (int c = 0; c < nsmall; c++) if (csz[c] != small[c]) return 0; return 1; }
    int x = el[i], opened = 0;
    for (int c = 0; c < nsmall; c++) {
        if (csz[c] == 0) { if (opened) continue; /* empty classes of equal size are interchangeable */ }
        if (csz[c] >= small[c]) continue;
        if (csz[c] == 0) { int dup = 0; for (int e = 0; e < c; e++) if (csz[e] == 0 && small[e] == small[c]) dup = 1; if (dup) continue; opened = 1; }
        u128 d = 0; int ok = 1;
        for (int j = 0; j < csz[c]; j++) { u128 b = (u128)1 << (x - cmem[c][j]); if ((cdif[c] | d) & b) { ok = 0; break; } d |= b; }
        if (!ok) continue;
        cmem[c][csz[c]++] = x; cdif[c] |= d;
        int r = sp(i + 1, exact);
        csz[c]--; cdif[c] &= ~d;
        if (r) return 1;
    }
    return 0;
}
static int splits(u128 m, int exact) {
    splitcalls++;
    ne = 0; for (int x = 1; x <= N; x++) if ((m >> x) & 1) el[ne++] = x;
    if (ne > smalltot) return 0;
    memset(csz, 0, sizeof csz); memset(cdif, 0, sizeof cdif);
    return sp(0, exact);
}

static unsigned long long rs = 88172645463325252ULL;
static unsigned long long xr(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return rs; }
static double est_total, est_prod; static int est_mode;

// returns via *bx the chosen number and fills cnt; handles forced-leftover
static void rec(int d, int *L, long nL, u128 unc, u128 left, int k10) {
    if (timeout) return;
    nodes++; nodes_at[d]++;
    if (now() > deadline) { timeout = 1; return; }
    if (k10 == 0) {  // all remaining numbers go to the leftover
        if (splits(left | unc, 1)) { sols++; if (MODE == 0) timeout = 2;
            if (!printed) { printed = 1; fprintf(stderr, "SOL");
                for (int c = 0; c < nch; c++) { fprintf(stderr, " ["); int f1 = 1; for (int x = 1; x <= N; x++) if ((chosen[c] >> x) & 1) { fprintf(stderr, f1 ? "%d" : ",%d", x); f1 = 0; } fprintf(stderr, "]"); }
                { u128 lf = left | unc; fprintf(stderr, " LEFT["); int f1 = 1; for (int x = 1; x <= N; x++) if ((lf >> x) & 1) { fprintf(stderr, f1 ? "%d" : ",%d", x); f1 = 0; } fprintf(stderr, "]"); }
                fprintf(stderr, "\n"); fflush(stderr); } }
        return;
    }
    int cnt[130] = {0};
    for (long i = 0; i < nL; i++) {
        unsigned long long lo = (unsigned long long)M[L[i]], hi = (unsigned long long)(M[L[i]] >> 64);
        while (lo) { cnt[__builtin_ctzll(lo)]++; lo &= lo - 1; }
        while (hi) { cnt[64 + __builtin_ctzll(hi)]++; hi &= hi - 1; }
    }
    // forced: uncovered numbers in no remaining 10-set must be leftover
    u128 forced = 0;
    for (int x = 1; x <= N; x++) if (((unc >> x) & 1) && cnt[x] == 0) forced |= (u128)1 << x;
    if (forced) {
        u128 nl = left | forced;
        if (pc(nl) > smalltot || !splits(nl, 0)) return;
        // forced numbers never appear in L, so L stays valid
        nodes_at[d]--; nodes--;
        rec(d, L, nL, unc & ~forced, nl, k10);
        return;
    }
    // remaining sets must cover all uncovered numbers except at most (smalltot - |left|)
    int best = -1, bc = 1 << 30, bleft = 0;
    for (int x = 1; x <= N; x++) if ((unc >> x) & 1) {
        int o = cnt[x];
        if (o >= bc) continue;                // cannot beat even with leftover option 0
        int lo = (pc(left) < smalltot) && splits(left | ((u128)1 << x), 0);
        if (o + lo < bc) { bc = o + lo; best = x; bleft = lo; }
    }
    if (bc == 0) return;
    u128 bx = (u128)1 << best;
    if (est_mode) {  // one random child
        long pick = xr() % bc;
        est_prod *= bc; est_total += est_prod;
        if (pick < cnt[best]) {
            long seen = 0;
            for (long i = 0; i < nL; i++) if (M[L[i]] & bx) { if (seen++ == pick) {
                u128 c = M[L[i]]; long k = 0; int *NL = lev[d + 1];
                for (long j = 0; j < nL; j++) if (!(M[L[j]] & c)) NL[k++] = L[j];
                rec(d + 1, NL, k, unc & ~c, left, k10 - 1); return; } }
        } else {
            long k = 0; int *NL = lev[d + 1];
            for (long j = 0; j < nL; j++) if (!(M[L[j]] & bx)) NL[k++] = L[j];
            rec(d + 1, NL, k, unc & ~bx, left | bx, k10);
        }
        return;
    }
    for (long i = 0; i < nL; i++) {
        u128 c = M[L[i]];
        if (!(c & bx)) continue;
        long k = 0; int *NL = lev[d + 1];
        for (long j = 0; j < nL; j++) if (!(M[L[j]] & c)) NL[k++] = L[j];
        chosen[nch++] = c;
        rec(d + 1, NL, k, unc & ~c, left, k10 - 1);
        nch--;
        if (timeout) return;
    }
    if (bleft) {
        long k = 0; int *NL = lev[d + 1];
        for (long j = 0; j < nL; j++) if (!(M[L[j]] & bx)) NL[k++] = L[j];
        rec(d + 1, NL, k, unc & ~bx, left | bx, k10);
    }
}


static u128 *ALL; static long nALL; static long levcap;
// Prefix cache: the fixed classes arrive as groups of 11 numbers, and consecutive cases
// share leading groups (tuples are in lexicographic order). Level g holds the sets
// disjoint from groups 0..g; it is rebuilt only from the first group that changed.
static u128 *PL[8]; static long PN[8]; static u128 PK[8]; static int PG = 0;
static void reset_and_solve_groups(u128 *grp, int ng, u128 unc, int k10, const char *sm, double secs, char *out, int reuse) {
    nsmall = 0; smalltot = 0;
    if (strcmp(sm, "-")) for (const char *p = sm; *p; ) { char *e; small[nsmall] = strtol(p, &e, 10); p = e; smalltot += small[nsmall++]; if (*p == ',') p++; }
    if (!reuse) {
        int g0 = 0; while (g0 < PG && g0 < ng && PK[g0] == grp[g0]) g0++;
        for (int g = g0; g < ng; g++) {
            u128 *src = g ? PL[g - 1] : ALL; long sn = g ? PN[g - 1] : nALL; long k = 0;
            for (long i = 0; i < sn; i++) if (!(src[i] & grp[g])) PL[g][k++] = src[i];
            PN[g] = k; PK[g] = grp[g];
        }
        PG = ng;
    }
    nM = 0;
    u128 *src = PL[ng - 1]; long sn = PN[ng - 1];
    for (long i = 0; i < sn; i++) if (!(src[i] & ~unc)) M[nM++] = src[i];
    if (nM + 1 > levcap) { levcap = 2 * (nM + 1); for (int d = 0; d < 64; d++) lev[d] = realloc(lev[d], levcap * sizeof(int)); }
    for (long i = 0; i < nM; i++) lev[0][i] = i;
    nch = 0; printed = 0;
    nodes = sols = splitcalls = 0; timeout = 0; memset(nodes_at, 0, sizeof nodes_at);
    double ts = now(); deadline = ts + secs;
    rec(0, lev[0], nM, unc, 0, k10);
    sprintf(out, "%s k10=%d small=%s sets=%ld sols=%lld nodes=%lld secs=%.2f",
            timeout == 1 ? "TIMEOUT" : (sols ? "FOUND" : "EXHAUSTED"), k10, sm, nM, sols, nodes, now() - ts);
}
// xcover3 N SETS10.bin batch PROFILES SECS < lines "NAME n1,n2,..."
// PROFILES like "4:8/3:9,9"; MODE fixed to "count all".
static int batch_main(int argc, char **argv) {
    double secs = atof(argv[5]); MODE = 1;
    FILE *f = fopen(argv[2], "rb"); unsigned char b[SZ]; long cap = 1 << 20; ALL = malloc(cap * sizeof *ALL);
    while (fread(b, 1, SZ, f) == SZ) {
        u128 m = 0; int ok = 1;
        for (int i = 0; i < SZ; i++) { if (b[i] > N) { ok = 0; break; } m |= (u128)1 << b[i]; }
        if (!ok) continue;
        if (nALL == cap) { cap *= 2; ALL = realloc(ALL, cap * sizeof *ALL); }
        ALL[nALL++] = m;
    }
    fclose(f);
    M = malloc((nALL + 1) * sizeof *M); fprintf(stderr, "loaded %ld sets of size %d\n", nALL, SZ);
    for (int g = 0; g < 8; g++) PL[g] = NULL;
    lev = calloc(64, sizeof *lev);
    char prof[256]; strcpy(prof, argv[4]);
    int pk[8]; char psm[8][32]; int np = 0;
    for (char *t = strtok(prof, "/"); t; t = strtok(0, "/")) { pk[np] = atoi(t); strcpy(psm[np], strchr(t, ':') + 1); np++; }
    char name[256], list[1024], out[512];
    while (scanf("%255s %1023s", name, list) == 2) {
        u128 fixed = 0, grp[8] = {0}; int cntx = 0;
        for (char *p = list; *p; ) { long x = strtol(p, &p, 10); fixed |= (u128)1 << x; grp[cntx / GS] |= (u128)1 << x; cntx++; if (*p == ',') p++; }
        int ng = (cntx + GS - 1) / GS;
        for (int g = 0; g < ng; g++) if (!PL[g]) PL[g] = malloc((nALL + 1) * sizeof(u128));
        u128 unc = 0; for (int i = 1; i <= N; i++) if (!((fixed >> i) & 1)) unc |= (u128)1 << i;
        printf("%s", name);
        for (int q = 0; q < np; q++) { reset_and_solve_groups(grp, ng, unc, pk[q], psm[q], secs, out, q > 0); printf(" | %s", out); }
        printf("\n"); fflush(stdout);
    }
    return 0;
}
int main(int argc, char **argv) {
    if (argc >= 6 && !strcmp(argv[3], "batch")) { N = atoi(argv[1]); return batch_main(argc, argv); }
    if (argc < 7) { fprintf(stderr, "uso: %s N SETS10.bin K10 SMALL MODE SECS < fixed\n", argv[0]); return 2; }
    N = atoi(argv[1]); K10 = atoi(argv[3]); MODE = atoi(argv[5]); double secs = atof(argv[6]);
    if (strcmp(argv[4], "-")) for (char *p = argv[4]; *p; ) { small[nsmall] = strtol(p, &p, 10); smalltot += small[nsmall++]; if (*p == ',') p++; }
    u128 fixed = 0; int x; while (scanf(" %d ,", &x) == 1) fixed |= (u128)1 << x;
    u128 unc = 0; int nfree = 0;
    for (int i = 1; i <= N; i++) if (!((fixed >> i) & 1)) { unc |= (u128)1 << i; nfree++; }
    if (nfree != K10 * SZ + smalltot) { fprintf(stderr, "free=%d != %d*10+%d\n", nfree, K10, smalltot); return 2; }
    FILE *f = fopen(argv[2], "rb"); unsigned char b[SZ]; long cap = 1 << 20; M = malloc(cap * sizeof *M);
    while (fread(b, 1, SZ, f) == SZ) {
        u128 m = 0; int ok = 1;
        for (int i = 0; i < SZ; i++) { if (b[i] > N) { ok = 0; break; } m |= (u128)1 << b[i]; }
        if (!ok || (m & ~unc)) continue;
        if (nM == cap) { cap *= 2; M = realloc(M, cap * sizeof *M); }
        M[nM++] = m;
    }
    fclose(f);
    lev = malloc(64 * sizeof *lev);
    for (int d = 0; d < 64; d++) lev[d] = malloc((nM + 1) * sizeof(int));
    for (long i = 0; i < nM; i++) lev[0][i] = i;
    t0 = now(); deadline = t0 + secs;
    if (MODE >= 2) {
        est_mode = 1; double sum = 0, sum2 = 0;
        for (int p = 0; p < MODE; p++) { est_total = 1; est_prod = 1; rec(0, lev[0], nM, unc, 0, K10); sum += est_total; sum2 += est_total * est_total; }
        double m = sum / MODE, sd = sqrt(fmax(0, sum2 / MODE - m * m) / MODE);
        printf("ESTIMATE k10=%d small=%s sets=%ld nodes=%.3g +- %.2g (probes=%d, %.1fs, %.3g ms/node)\n", K10, argv[4], nM, m, sd, MODE, now() - t0, 1000 * (now() - t0) / nodes);
        return 0;
    }
    rec(0, lev[0], nM, unc, 0, K10);
    printf("%s k10=%d small=%s sets=%ld sols=%lld nodes=%lld splits=%lld secs=%.1f\n",
           timeout == 1 ? "TIMEOUT" : (sols ? "FOUND" : "EXHAUSTED"), K10, argv[4], nM, sols, nodes, splitcalls, now() - t0);
    return 0;
}
