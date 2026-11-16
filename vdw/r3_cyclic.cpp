// exact r_3(Z_q): largest subset of Z_q with no x, x+d, x+2d (d != 0 mod q) all in it; q <= 64
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <vector>
int q, best; std::vector<std::vector<std::pair<int,int>>> pairs; // for element z: pairs (a,b) with {a,b,z} an AP
void dfs(int i, uint64_t in, uint64_t banned, int size) {
    if (size + (q - i) <= best) return;
    if (i == q) { best = size; return; }
    if (!(banned >> i & 1)) {                     // take i
        uint64_t nb = banned;
        bool ok = true;
        for (auto [a, b] : pairs[i]) {            // with i in, a in -> ban b, etc.
            bool ina = in >> a & 1, inb = in >> b & 1;
            if (a == b) { if (ina) { ok = false; break; } else nb |= 1ULL << a; continue; }
            if (ina && inb) { ok = false; break; }
            if (ina) nb |= 1ULL << b;
            if (inb) nb |= 1ULL << a;
        }
        if (ok) dfs(i + 1, in | 1ULL << i, nb, size + 1);
    }
    dfs(i + 1, in, banned, size);
}
int main(int argc, char** argv) {
    int lo = atoi(argv[1]), hi = atoi(argv[2]);
    for (q = lo; q <= hi; ++q) {
        pairs.assign(q, {});
        for (int x = 0; x < q; ++x) for (int d = 1; d < q; ++d) {
            int y = (x + d) % q, z = (x + 2 * d) % q;
            int t[3] = {x, y, z};
            for (int k = 0; k < 3; ++k) pairs[t[k]].push_back({t[(k + 1) % 3], t[(k + 2) % 3]});
        }
        best = 0; dfs(0, 0, 0, 0);
        printf("%d %d\n", q, best); fflush(stdout);
    }
}
