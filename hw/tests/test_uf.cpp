// Small self-contained checks for the union-find decoder (no framework).
#include "../src/uf.hpp"
#include <cstdio>

static int fails = 0;
#define CHECK(c) do { if (!(c)) { std::printf("FAIL %s:%d %s\n", __FILE__, __LINE__, #c); ++fails; } } while (0)

// repetition code: detectors 0..n-1 on a line, boundary at both ends.
static uf::Graph line(int n, int obs_edge) {
    uf::Graph g; g.num_det = n; g.num_edges = n + 1;
    for (int e = 0; e <= n; ++e) {
        int u = e == 0 ? 0 : e - 1, v = e == 0 ? n : (e == n ? n : e);
        if (e == n) u = n - 1;
        g.eu.push_back(u); g.ev.push_back(v); g.eobs.push_back(e == obs_edge);
    }
    g.build_adjacency();
    return g;
}

int main() {
    uf::Graph g = line(5, 0);  // edge 0: boundary-det0 carries the observable
    uf::Decoder d(g);
    uint8_t s0[5] = {0,0,0,0,0};
    CHECK(d.decode(s0, nullptr) == 0);
    uint8_t s1[5] = {1,0,0,0,0};          // single defect near left boundary -> left edge, flips obs
    CHECK(d.decode(s1, nullptr) == 1);
    uint8_t s2[5] = {0,0,0,0,1};          // near right boundary -> right edge, no obs
    CHECK(d.decode(s2, nullptr) == 0);
    uint8_t s3[5] = {0,1,1,0,0};          // pair in the middle
    CHECK(d.decode(s3, nullptr) == 0);
    // erasure: defect at det 2, erase edges 0,1,2 (boundary..det2) -> path through obs edge
    uint8_t s4[5] = {0,0,1,0,0};
    uint8_t er[6] = {1,1,1,0,0,0};
    CHECK(d.decode(s4, er) == 1);
    uint8_t er2[6] = {0,0,0,1,1,1};       // erase right side instead -> no obs flip
    CHECK(d.decode(s4, er2) == 0);
    // decoder state is reset between shots
    CHECK(d.decode(s0, nullptr) == 0);
    std::printf(fails ? "uf tests FAILED\n" : "uf tests passed\n");
    return fails ? 1 : 0;
}
