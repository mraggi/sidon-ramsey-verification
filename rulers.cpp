// rulers M L T OUT : write every M-mark Golomb ruler 0 = a_0 < ... < a_{M-1} <= L (both
// orientations; M bytes per ruler, sorted) to OUT and print the number found per length.
//
// T = "t1,t2,...,t_{M-1}" are lower bounds for the length of an r-mark Golomb ruler, r < M
// (t1 = 0, t2 = 1). They are only used to cut branches, so the output is complete as long as
// each t_r is at most the true optimal length g_r. Running this program for r marks with
// L = x - 1 and getting no ruler proves g_r >= x, given the bounds for fewer marks.
//
// Search: the marks are placed in increasing order. Three bit sets of distances describe the
// marks placed so far, with p the last one: D, the differences used; P, the distances from p
// back to the earlier marks; and C, the gaps h for which a mark at p + h would repeat a
// difference, that is, h + e is in D for e = 0 or some e in P. After a new mark at p + g,
// P' = (P << g) | {g}, D' = D | P' and C' = (C >> g) | D'. (A gap h with h + e in P' but not
// in C >> g or D' cannot exist: h + e = p' - a_j and e = p' - a_i give h = a_i - a_j in D'.)
// Two cuts, both valid for every Golomb ruler: (1) the marks a_j..a_{M-1} form a ruler with
// M - j marks, so a_{M-1} - a_j >= t_{M-j}; (2) the M-1-j gaps after a_j are distinct
// differences not used yet, so their sum a_{M-1} - a_j is at least the sum of the M-1-j
// smallest unused differences. Each ruler is found once per mirror pair: the mirror image
// a'_i = l - a_{M-1-i} of a ruler of length l has middle sum s' = 2l - s, where
// s = a_lo + a_hi with lo = floor((M-1)/2), hi = floor(M/2); the search keeps the image with
// s < l (or, if s = l, the lexicographically smaller one) and writes both. Since s <= l <= L,
// branches with s > L are cut as soon as a_hi is placed.
//
//   g++ -O3 -march=native -std=c++17 -fopenmp -o rulers rulers.cpp
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

typedef uint64_t u64;
static int M, L, T[64], LO, HI;

struct Bits {
    u64 w[4] = {0, 0, 0, 0};
    bool has(int d) const { return (w[d >> 6] >> (d & 63)) & 1; }
    void add(int d) { w[d >> 6] |= u64(1) << (d & 63); }
};
static inline Bits shl(const Bits &b, int g) {  // 0 < g < 64
    Bits r;
    r.w[3] = (b.w[3] << g) | (b.w[2] >> (64 - g));
    r.w[2] = (b.w[2] << g) | (b.w[1] >> (64 - g));
    r.w[1] = (b.w[1] << g) | (b.w[0] >> (64 - g));
    r.w[0] = b.w[0] << g;
    return r;
}
static inline Bits shr(const Bits &b, int g) {  // 0 < g < 64
    Bits r;
    r.w[0] = (b.w[0] >> g) | (b.w[1] << (64 - g));
    r.w[1] = (b.w[1] >> g) | (b.w[2] << (64 - g));
    r.w[2] = (b.w[2] >> g) | (b.w[3] << (64 - g));
    r.w[3] = b.w[3] >> g;
    return r;
}
static inline Bits shift_left(Bits b, int g) { while (g >= 63) { b = shl(b, 63); g -= 63; } return g ? shl(b, g) : b; }
static inline Bits shift_right(Bits b, int g) { while (g >= 63) { b = shr(b, 63); g -= 63; } return g ? shr(b, g) : b; }
static inline Bits orb(const Bits &a, const Bits &b) { Bits r; for (int k = 0; k < 4; k++) r.w[k] = a.w[k] | b.w[k]; return r; }

// sum of the r smallest positive integers not in b, stopping once it exceeds cap
static inline int smallest_unused_sum(const Bits &b, int r, int cap) {
    int s = 0;
    for (int k = 0; k < 4 && r > 0; k++) {
        u64 free = ~b.w[k];
        if (k == 0) free &= ~u64(1);   // 0 is not a difference
        while (free && r > 0) {
            s += k * 64 + __builtin_ctzll(free);
            if (s > cap) return s;
            free &= free - 1;
            r--;
        }
    }
    return r > 0 ? cap + 1 : s;
}

static void emit(const int *a, std::vector<uint8_t> &out) {
    int l = a[M - 1], s = a[LO] + a[HI];
    if (s > l) return;                       // the mirror image is the one kept
    uint8_t r[64], m[64];
    for (int i = 0; i < M; i++) { r[i] = a[i]; m[i] = l - a[M - 1 - i]; }
    if (s == l && memcmp(m, r, M) < 0) return;  // both images have s = l: keep the smaller
    out.insert(out.end(), r, r + M);
    out.insert(out.end(), m, m + M);
}

static void dfs(int j, int *a, const Bits &P, const Bits &D, const Bits &C, std::vector<uint8_t> &out) {
    // marks a_0..a_{j-1} placed; place a_j = p + g
    int p = a[j - 1], after = M - 1 - j;
    int gmax = L - p - T[after + 1];         // cut (1)
    for (int g = 1; g <= gmax; g++) {
        if (C.has(g)) continue;
        Bits nP = shift_left(P, g); nP.add(g);
        Bits nD = orb(D, nP);
        int q = p + g;
        if (after > 0 && smallest_unused_sum(nD, after, L - q) > L - q) continue;  // cut (2)
        a[j] = q;
        if (j == HI && a[LO] + a[HI] > L) continue;  // mirror: s > L >= l
        if (after == 0) { emit(a, out); continue; }
        Bits nC = orb(shift_right(C, g), nD);
        dfs(j + 1, a, nP, nD, nC, out);
    }
}

int main(int argc, char **argv) {
    if (argc != 5) { fprintf(stderr, "usage: %s M L t1,t2,...,t_{M-1} OUT\n", argv[0]); return 2; }
    M = atoi(argv[1]); L = atoi(argv[2]);
    if (M < 3 || M > 40 || L < 1 || L > 200) { fprintf(stderr, "need 3 <= M <= 40, L <= 200\n"); return 2; }
    std::vector<int> t;
    for (char *s = strtok(argv[3], ","); s; s = strtok(nullptr, ",")) t.push_back(atoi(s));
    if ((int)t.size() != M - 1 || t[0] != 0 || t[1] != 1) { fprintf(stderr, "need M-1 bounds, t1=0, t2=1\n"); return 2; }
    for (int r = 1; r < M; r++) T[r] = t[r - 1];
    LO = (M - 1) / 2; HI = M / 2;
    // one task per first gap a_1 and second gap
    std::vector<std::pair<int, int>> tasks;
    for (int a1 = 1; a1 <= L; a1++)
        for (int g2 = 1; a1 + g2 <= L; g2++)
            if (g2 != a1) tasks.push_back({a1, g2});
    std::vector<std::vector<uint8_t>> found(tasks.size());
    #pragma omp parallel for schedule(dynamic, 1)
    for (long k = 0; k < (long)tasks.size(); k++) {
        int a[64] = {0};
        Bits P, D, C;
        // place a_1 = g1 and a_2 = a_1 + g2 through the same recurrence as dfs, so every check applies
        int g1 = tasks[k].first, g2 = tasks[k].second;
        a[1] = g1;
        if (g1 > L - T[M - 1]) continue;
        P.add(g1); D = P; C = D;                       // after marks {0, g1}
        if (M > 3 && smallest_unused_sum(D, M - 2, L - g1) > L - g1) continue;
        if (HI == 1 && a[LO] + a[HI] > L) continue;
        if (C.has(g2) || a[1] + g2 > L - T[M - 2]) continue;
        Bits nP = shift_left(P, g2); nP.add(g2);
        Bits nD = orb(D, nP);
        a[2] = g1 + g2;
        if (M > 3 && smallest_unused_sum(nD, M - 3, L - a[2]) > L - a[2]) continue;
        if (HI == 2 && a[LO] + a[HI] > L) continue;
        if (M == 3) { emit(a, found[k]); continue; }
        Bits nC = orb(shift_right(C, g2), nD);
        dfs(3, a, nP, nD, nC, found[k]);
    }
    std::vector<std::vector<uint8_t>> all;
    for (auto &f : found)
        for (size_t i = 0; i < f.size(); i += M) all.emplace_back(f.begin() + i, f.begin() + i + M);
    std::sort(all.begin(), all.end());
    all.erase(std::unique(all.begin(), all.end()), all.end());
    std::vector<long> per_len(L + 1, 0);
    FILE *fo = fopen(argv[4], "wb");
    for (auto &r : all) { fwrite(r.data(), 1, M, fo); per_len[r[M - 1]]++; }
    fclose(fo);
    printf("{\"M\": %d, \"L\": %d, \"count\": %zu, \"by_length\": {", M, L, all.size());
    bool first = true;
    for (int l = 0; l <= L; l++) if (per_len[l]) { printf("%s\"%d\": %ld", first ? "" : ", ", l, per_len[l]); first = false; }
    printf("}}\n");
    return 0;
}
