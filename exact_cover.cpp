// exact_cover U.txt SETS.bin M F : count the ways to write the point set U (one line, points separated
// by spaces) as a disjoint union of members of the list SETS.bin (M bytes per set, all inside U)
// plus F further points that form a Sidon set. Prints the number of solutions and the first one.
//
// Exhaustive search (Knuth's Algorithm X with the leftover points as an extra option): take the
// undecided point with the fewest options, where its options are the live sets containing it
// (disjoint from everything decided so far) and, if fewer than F points are left over and the
// leftover set stays Sidon with it, leaving it over. Every solution is found exactly once.
//
// Used by verify_upper.py --exhaustive, as a check of SR(10) <= 132 that does not use the weights.
//
//   g++ -O3 -march=native -std=c++17 -o exact_cover exact_cover.cpp
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

typedef uint64_t u64;
static int n, M, F, W;                       // n points of U, W words per candidate bitset
static std::vector<int> pts;                 // the points of U
static std::vector<std::vector<int>> sets;   // candidate sets, as indices into pts
static std::vector<std::vector<u64>> contains;  // contains[x] = bitset of candidates containing point x
static long solutions = 0, nodes = 0;
static std::vector<int> chosen, first_solution_sets, first_left, leftover;
static bool decided[256];

static bool leftover_sidon_with(int x) {     // is leftover + pts[x] still a Sidon set? (differences)
    std::vector<int> d;
    for (int i = 0; i < (int)leftover.size(); i++)
        for (int j = i + 1; j < (int)leftover.size(); j++) d.push_back(abs(pts[leftover[i]] - pts[leftover[j]]));
    for (int y : leftover) {
        int e = abs(pts[x] - pts[y]);
        for (int f : d) if (f == e) return false;
        d.push_back(e);
    }
    return true;
}

static void search(std::vector<u64> &live) {
    nodes++;
    int best = -1, bestc = 1 << 30; bool best_left = false;
    for (int x = 0; x < n; x++) {
        if (decided[x]) continue;
        int c = 0;
        for (int k = 0; k < W; k++) c += __builtin_popcountll(live[k] & contains[x][k]);
        bool can_left = (int)leftover.size() < F && leftover_sidon_with(x);
        c += can_left;
        if (c < bestc) { bestc = c; best = x; best_left = can_left; if (c == 0) break; }
    }
    if (best < 0) {                           // every point decided
        if ((int)leftover.size() == F) {
            if (solutions == 0) { first_solution_sets = chosen; first_left = leftover; }
            solutions++;
        }
        return;
    }
    if (bestc == 0) return;
    int x = best;
    if (best_left) {                          // leave x over
        decided[x] = true; leftover.push_back(x);
        std::vector<u64> nl(live);
        for (int k = 0; k < W; k++) nl[k] &= ~contains[x][k];
        search(nl);
        leftover.pop_back(); decided[x] = false;
    }
    for (int k = 0; k < W; k++) {             // cover x by a live set
        u64 b = live[k] & contains[x][k];
        while (b) {
            int c = k * 64 + __builtin_ctzll(b); b &= b - 1;
            std::vector<u64> nl(live);
            for (int y : sets[c]) { decided[y] = true; for (int j = 0; j < W; j++) nl[j] &= ~contains[y][j]; }
            chosen.push_back(c);
            search(nl);
            chosen.pop_back();
            for (int y : sets[c]) decided[y] = false;
        }
    }
}

int main(int argc, char **argv) {
    if (argc != 5) { fprintf(stderr, "usage: %s U.txt SETS.bin M F\n", argv[0]); return 2; }
    M = atoi(argv[3]); F = atoi(argv[4]);
    std::ifstream in(argv[1]); std::string line; std::getline(in, line);
    std::istringstream ss(line); int p; int index[256]; for (int i = 0; i < 256; i++) index[i] = -1;
    while (ss >> p) { index[p] = pts.size(); pts.push_back(p); }
    n = pts.size();
    FILE *f = fopen(argv[2], "rb"); if (!f) { perror(argv[2]); return 2; }
    std::vector<uint8_t> b(M);
    while (fread(b.data(), 1, M, f) == (size_t)M) {
        std::vector<int> s;
        for (int i = 0; i < M; i++) { if (index[b[i]] < 0) { fprintf(stderr, "set not inside U\n"); return 2; } s.push_back(index[b[i]]); }
        sets.push_back(s);
    }
    fclose(f);
    W = (sets.size() + 63) / 64; if (W == 0) W = 1;
    contains.assign(n, std::vector<u64>(W, 0));
    for (size_t c = 0; c < sets.size(); c++) for (int y : sets[c]) contains[y][c >> 6] |= u64(1) << (c & 63);
    std::vector<u64> live(W, 0);
    for (size_t c = 0; c < sets.size(); c++) live[c >> 6] |= u64(1) << (c & 63);
    search(live);
    printf("points %d candidates %zu leftover %d solutions %ld nodes %ld\n", n, sets.size(), F, solutions, nodes);
    if (solutions) {
        for (int c : first_solution_sets) { printf("SET"); for (int y : sets[c]) printf(" %d", pts[y]); printf("\n"); }
        if (F) { printf("LEFT"); for (int y : first_left) printf(" %d", pts[y]); printf("\n"); }
    }
    return 0;
}
