// sets_inside RULERS.bin M N CASES.txt [OUTDIR] : for each line of CASES.txt, print the number
// of M-element Sidon subsets of [N] that lie inside the set U of the line, and the smallest
// weight of such a set; with OUTDIR, also write these sets to OUTDIR/<line number>.bin
// (M bytes per set, in increasing order, sorted).
// A line has N entries, one per point 1..N: an integer weight for a point of U and "x" for a
// point outside U. The M-element Sidon subsets of [N] are the
// translates a + s (s = 1, ..., N - length) of the M-mark Golomb rulers of length at most
// N - 1, which RULERS.bin lists (M bytes per ruler, first mark 0), as written by rulers.cpp.
//
//   g++ -O3 -march=native -std=c++17 -fopenmp -o sets_inside sets_inside.cpp
#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

typedef uint64_t u64;

int main(int argc, char **argv) {
    if (argc != 5 && argc != 6) { fprintf(stderr, "usage: %s RULERS.bin M N CASES.txt [OUTDIR]\n", argv[0]); return 2; }
    bool keep = argc == 6;
    int M = atoi(argv[2]), N = atoi(argv[3]);
    if (M < 1 || N < 1 || N > 255) { fprintf(stderr, "need N <= 255\n"); return 2; }
    std::vector<uint8_t> R;
    {
        FILE *f = fopen(argv[1], "rb");
        if (!f) { perror(argv[1]); return 2; }
        int c;
        while ((c = fgetc(f)) != EOF) R.push_back((uint8_t)c);
        fclose(f);
    }
    if (R.size() % M) { fprintf(stderr, "file size is not a multiple of M\n"); return 2; }
    long nR = R.size() / M;
    for (long r = 0; r < nR; r++)
        if (R[r * M] != 0 || R[r * M + M - 1] > N - 1) { fprintf(stderr, "ruler %ld not normalized or too long\n", r); return 2; }

    std::vector<std::vector<u64>> inU;   // bit masks of the sets U
    std::vector<std::vector<long>> W;     // weights, indexed by point
    std::ifstream in(argv[4]);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::vector<u64> m(4, 0);
        std::vector<long> w(N + 1, 0);
        std::string tok;
        for (int x = 1; x <= N; x++) {
            if (!(ss >> tok)) { fprintf(stderr, "line with fewer than N entries\n"); return 2; }
            if (tok == "x") continue;
            m[x >> 6] |= u64(1) << (x & 63);
            w[x] = atol(tok.c_str());
        }
        inU.push_back(m);
        W.push_back(w);
    }
    int nC = inU.size();
    std::vector<long> count(nC, 0), minw(nC, LONG_MAX);
    std::vector<std::vector<std::vector<uint8_t>>> kept(nC);
    #pragma omp parallel
    {
        std::vector<long> cnt(nC, 0), mn(nC, LONG_MAX);
        std::vector<std::vector<std::vector<uint8_t>>> kp(nC);
        #pragma omp for schedule(dynamic, 4096)
        for (long r = 0; r < nR; r++) {
            const uint8_t *a = &R[r * M];
            for (int s = 1; s + a[M - 1] <= N; s++) {
                u64 m[4] = {0, 0, 0, 0};
                for (int i = 0; i < M; i++) { int x = a[i] + s; m[x >> 6] |= u64(1) << (x & 63); }
                for (int c = 0; c < nC; c++) {
                    if ((m[0] & ~inU[c][0]) | (m[1] & ~inU[c][1]) | (m[2] & ~inU[c][2]) | (m[3] & ~inU[c][3])) continue;
                    long w = 0;
                    for (int i = 0; i < M; i++) w += W[c][a[i] + s];
                    cnt[c]++;
                    if (w < mn[c]) mn[c] = w;
                    if (keep) { std::vector<uint8_t> v(M); for (int i = 0; i < M; i++) v[i] = a[i] + s; kp[c].push_back(v); }
                }
            }
        }
        #pragma omp critical
        for (int c = 0; c < nC; c++) {
            count[c] += cnt[c];
            if (mn[c] < minw[c]) minw[c] = mn[c];
            kept[c].insert(kept[c].end(), kp[c].begin(), kp[c].end());
        }
    }
    if (keep)
        for (int c = 0; c < nC; c++) {
            std::sort(kept[c].begin(), kept[c].end());
            std::string fn = std::string(argv[5]) + "/" + std::to_string(c) + ".bin";
            FILE *f = fopen(fn.c_str(), "wb");
            if (!f) { perror(fn.c_str()); return 2; }
            for (auto &v : kept[c]) fwrite(v.data(), 1, M, f);
            fclose(f);
        }
    for (int c = 0; c < nC; c++) printf("%ld %ld\n", count[c], count[c] ? minw[c] : 0);
    return 0;
}
