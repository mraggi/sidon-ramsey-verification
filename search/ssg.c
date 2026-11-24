// all K-element Sidon sets in [N], pruned by optimal Golomb ruler lengths; binary output
// (K bytes per set, sorted) to ssK_N.bin.  OpenMP over the first two marks.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned __int128 u128;
static const int G[] = {0,0,1,3,6,11,17,25,34,44,55,72,85,106,127};
static int N, K;
typedef struct { unsigned char *buf; long n, cap; } Out;
static void rec(int k, int *a, u128 diffs, Out *o) {
    if (k == K) { if (o->n == o->cap) o->buf = realloc(o->buf, (o->cap = o->cap ? 2*o->cap : 4096) * K);
        for (int i = 0; i < K; i++) o->buf[o->n*K+i] = a[i]; o->n++; return; }
    int last = a[k-1], lim = N - G[K-k+1] + (k ? 0 : 0);
    // remaining K-k marks after 'last': ruler last..end has K-k+1 marks
    for (int x = last + 1; x + G[K-k] <= N; x++) {
        u128 d = 0; int ok = 1;
        for (int i = 0; i < k; i++) { int e = x - a[i]; u128 b = (u128)1 << e; if ((diffs | d) & b) { ok = 0; break; } d |= b; }
        if (ok) { a[k] = x; rec(k + 1, a, diffs | d, o); }
    }
    (void)lim;
}
int main(int argc, char **argv) {
    N = atoi(argv[1]); K = atoi(argv[2]);
    int pairs[20000][2], np = 0;
    for (int x = 1; x + G[K] <= N; x++) for (int y = x + 1; y + G[K-1] <= N; y++) { pairs[np][0] = x; pairs[np][1] = y; np++; }
    Out *outs = calloc(np, sizeof *outs);
    #pragma omp parallel for schedule(dynamic,1)
    for (int p = 0; p < np; p++) { int a[16]; a[0] = pairs[p][0]; a[1] = pairs[p][1];
        rec(2, a, (u128)1 << (a[1]-a[0]), &outs[p]); }
    char fn[64]; sprintf(fn, "ss%d_%d.bin", K, N); FILE *f = fopen(fn, "wb"); long tot = 0;
    for (int p = 0; p < np; p++) { fwrite(outs[p].buf, K, outs[p].n, f); tot += outs[p].n; }
    fclose(f); printf("%d-SS en [%d]: %ld\n", K, N, tot);
}
