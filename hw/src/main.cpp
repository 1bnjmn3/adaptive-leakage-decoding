// uf_decode graph.txt syndromes.u8 predictions.u8 [--erasures erasures.u8] [--repeat N]
// See hw/INTERFACE.md. Prints one JSON line with timing.
#include "uf.hpp"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

static std::vector<uint8_t> read_all(const char* path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) { std::fprintf(stderr, "cannot open %s\n", path); std::exit(2); }
    return std::vector<uint8_t>((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "usage: %s graph.txt syndromes.u8 predictions.u8 [--erasures f] [--repeat N]\n", argv[0]);
        return 2;
    }
    const char* er_path = nullptr;
    int repeat = 1;
    for (int i = 4; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--erasures") && i + 1 < argc) er_path = argv[++i];
        else if (!std::strcmp(argv[i], "--repeat") && i + 1 < argc) repeat = std::atoi(argv[++i]);
        else { std::fprintf(stderr, "bad arg %s\n", argv[i]); return 2; }
    }
    if (repeat < 1) repeat = 1;

    uf::Graph g;
    std::ifstream gf(argv[1]);
    if (!gf) { std::fprintf(stderr, "cannot open %s\n", argv[1]); return 2; }
    gf >> g.num_det >> g.num_edges;
    g.eu.resize(g.num_edges); g.ev.resize(g.num_edges); g.eobs.resize(g.num_edges);
    for (int32_t e = 0; e < g.num_edges; ++e) {
        long u, v, w, o;
        if (!(gf >> u >> v >> w >> o)) { std::fprintf(stderr, "bad graph line %d\n", e); return 2; }
        (void)w;  // unweighted union-find
        g.eu[e] = (int32_t)u;
        g.ev[e] = v < 0 ? g.num_det : (int32_t)v;
        g.eobs[e] = (uint8_t)(o & 1);
    }
    g.build_adjacency();

    const std::vector<uint8_t> syn = read_all(argv[2]);
    const size_t shots = g.num_det ? syn.size() / g.num_det : 0;
    std::vector<uint8_t> er;
    if (er_path) {
        er = read_all(er_path);
        if (er.size() != shots * (size_t)g.num_edges) { std::fprintf(stderr, "erasure size mismatch\n"); return 2; }
    }

    uf::Decoder dec(g);
    std::vector<uint8_t> pred(shots, 0);
    std::vector<int64_t> best(shots, INT64_MAX);
    for (int r = 0; r < repeat; ++r) {
        for (size_t s = 0; s < shots; ++s) {
            const uint8_t* sy = syn.data() + s * g.num_det;
            const uint8_t* e = er_path ? er.data() + s * g.num_edges : nullptr;
            const auto t0 = std::chrono::steady_clock::now();
            pred[s] = dec.decode(sy, e);
            const auto t1 = std::chrono::steady_clock::now();
            const int64_t ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
            if (ns < best[s]) best[s] = ns;
        }
    }
    int64_t total = 0, mx = 0;
    for (size_t s = 0; s < shots; ++s) { total += best[s]; if (best[s] > mx) mx = best[s]; }

    std::ofstream pf(argv[3], std::ios::binary);
    pf.write(reinterpret_cast<const char*>(pred.data()), (std::streamsize)pred.size());
    std::printf("{\"shots\":%zu,\"num_detectors\":%d,\"num_edges\":%d,\"total_decode_ns\":%lld,"
                "\"ns_per_shot\":%.3f,\"max_ns_per_shot\":%lld,\"static_bytes\":%zu}\n",
                shots, g.num_det, g.num_edges, (long long)total,
                shots ? (double)total / shots : 0.0, (long long)mx, dec.static_bytes());
    return 0;
}
