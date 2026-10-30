// Integer-only union-find decoder (Delfosse & Nickerson 2017) with erasure
// support (Delfosse & Zemor 2020): heralded edges start fully grown.
//
// Unweighted: every edge needs two half-edge growth steps. All memory is
// allocated once in init() from the graph size; decode() does no allocation.
#pragma once
#include <cstdint>
#include <cstddef>
#include <vector>

namespace uf {

struct Graph {
    int32_t num_det = 0;       // detectors 0..num_det-1; node num_det = boundary
    int32_t num_edges = 0;
    std::vector<int32_t> eu, ev;   // endpoints (boundary mapped to num_det)
    std::vector<uint8_t> eobs;     // flips observable 0?
    // CSR adjacency over nodes 0..num_det (boundary included)
    std::vector<int32_t> adj_off, adj_edge;
    void build_adjacency();
};

class Decoder {
public:
    explicit Decoder(const Graph& g);
    // syndrome: num_det bytes (0/1). erased: num_edges bytes or nullptr.
    // Returns predicted flip of observable 0.
    uint8_t decode(const uint8_t* syndrome, const uint8_t* erased);
    size_t static_bytes() const { return bytes_; }

private:
    const Graph& g_;
    int32_t n_;                 // nodes incl. boundary
    int32_t B_;                 // boundary node index
    // per node
    std::vector<int32_t> parent_, size_;
    std::vector<uint8_t> parity_, touches_bnd_, defect_, visited_;
    std::vector<int32_t> queue_, tree_edge_;
    // per edge
    std::vector<uint8_t> growth_;   // 0,1,2 half-edges grown
    std::vector<int32_t> fused_;    // edges that became fully grown this step
    std::vector<int32_t> grown_list_;  // edges with growth > 0 (for O(active) reset)
    // cluster member lists (singly linked, O(1) concatenation on union)
    std::vector<int32_t> next_, tail_;
    std::vector<int32_t> odd_roots_, odd_tmp_;  // active (odd, non-boundary) clusters
    std::vector<int32_t> dirty_;               // nodes whose state changed this shot
    std::vector<uint8_t> is_dirty_;
    int32_t n_grown_ = 0, n_dirty_ = 0;
    void mark(int32_t x) { if (!is_dirty_[x]) { is_dirty_[x] = 1; dirty_[n_dirty_++] = x; } }
    void grow(int32_t e, uint8_t inc, int32_t& nf);
    size_t bytes_ = 0;

    int32_t find(int32_t x);
    void unite(int32_t a, int32_t b);
    bool odd(int32_t root) const { return parity_[root] && !touches_bnd_[root]; }
};

}  // namespace uf
