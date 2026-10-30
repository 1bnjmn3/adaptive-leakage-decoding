#include "uf.hpp"

namespace uf {

void Graph::build_adjacency() {
    const int32_t n = num_det + 1;
    adj_off.assign(n + 1, 0);
    for (int32_t e = 0; e < num_edges; ++e) { adj_off[eu[e] + 1]++; adj_off[ev[e] + 1]++; }
    for (int32_t i = 0; i < n; ++i) adj_off[i + 1] += adj_off[i];
    adj_edge.assign(2 * num_edges, 0);
    std::vector<int32_t> fill(adj_off.begin(), adj_off.end() - 1);
    for (int32_t e = 0; e < num_edges; ++e) { adj_edge[fill[eu[e]]++] = e; adj_edge[fill[ev[e]]++] = e; }
}

Decoder::Decoder(const Graph& g) : g_(g), n_(g.num_det + 1), B_(g.num_det) {
    parent_.assign(n_, 0); size_.assign(n_, 1);
    parity_.assign(n_, 0); touches_bnd_.assign(n_, 0); defect_.assign(n_, 0); visited_.assign(n_, 0);
    queue_.assign(n_, 0); tree_edge_.assign(n_, -1);
    growth_.assign(g.num_edges, 0); fused_.assign(2 * g.num_edges, 0); grown_list_.assign(g.num_edges, 0);
    next_.assign(n_, -1); tail_.assign(n_, 0);
    odd_roots_.assign(n_, 0); odd_tmp_.assign(n_, 0);
    dirty_.assign(n_, 0); is_dirty_.assign(n_, 0);
    for (int32_t i = 0; i < n_; ++i) { parent_[i] = i; tail_[i] = i; }
    touches_bnd_[B_] = 1;
    // Everything above is allocated once here; decode() never allocates.
    // Count: per node 9 int32 arrays (parent,size,queue,tree_edge,next,tail,odd_roots,odd_tmp,dirty) + 5 flag bytes; per edge growth(1)
    // + fused(2 int32) + grown_list(int32); plus read-only graph tables.
    bytes_ = (size_t)n_ * (9 * sizeof(int32_t) + 5)
           + (size_t)g.num_edges * (1 + 3 * sizeof(int32_t))
           + (size_t)g.num_edges * (2 * sizeof(int32_t) + 1)             // eu, ev, eobs
           + (size_t)(n_ + 1 + 2 * g.num_edges) * sizeof(int32_t);       // CSR adjacency
}

int32_t Decoder::find(int32_t x) {
    int32_t r = x;
    while (parent_[r] != r) r = parent_[r];
    while (parent_[x] != r) { int32_t nx = parent_[x]; parent_[x] = r; x = nx; }  // path compression
    return r;
}

void Decoder::unite(int32_t a, int32_t b) {
    a = find(a); b = find(b);
    if (a == b) return;
    if (size_[a] < size_[b]) { int32_t t = a; a = b; b = t; }  // union by size
    mark(a); mark(b);
    parent_[b] = a;
    size_[a] += size_[b];
    parity_[a] ^= parity_[b];
    touches_bnd_[a] |= touches_bnd_[b];
    next_[tail_[a]] = b;      // append b's member list to a's
    tail_[a] = tail_[b];
}

void Decoder::grow(int32_t e, uint8_t inc, int32_t& nf) {
    if (growth_[e] == 0) grown_list_[n_grown_++] = e;
    const int32_t gnew = growth_[e] + inc;
    growth_[e] = (uint8_t)(gnew > 2 ? 2 : gnew);
    if (growth_[e] == 2) fused_[nf++] = e;
}

uint8_t Decoder::decode(const uint8_t* syn, const uint8_t* erased) {
    const int32_t E = g_.num_edges;
    const int32_t* eu = g_.eu.data();
    const int32_t* ev = g_.ev.data();
    const int32_t* off = g_.adj_off.data();
    const int32_t* adj = g_.adj_edge.data();
    n_grown_ = 0; n_dirty_ = 0;
    // 1. defects
    int32_t n_odd = 0;
    for (int32_t i = 0; i < g_.num_det; ++i) {
        if (syn[i]) { defect_[i] = 1; parity_[i] = 1; mark(i); odd_roots_[n_odd++] = i; }
    }
    // 2. erasures: heralded edges start fully grown (Delfosse-Zemor)
    if (erased) {
        int32_t nf = 0;
        for (int32_t e = 0; e < E; ++e)
            if (erased[e]) { growth_[e] = 0; grow(e, 2, nf); unite(eu[e], ev[e]); }
    }
    // 3. growth: every odd cluster grows half an edge along its whole boundary per step
    while (true) {
        int32_t m = 0;
        for (int32_t k = 0; k < n_odd; ++k) {          // refresh the active-root list
            const int32_t r = find(odd_roots_[k]);
            if (!odd(r) || is_dirty_[r] == 2) continue;  // 2 = already listed this step
            is_dirty_[r] = 2; odd_tmp_[m++] = r;
        }
        for (int32_t k = 0; k < m; ++k) is_dirty_[odd_tmp_[k]] = 1;
        if (m == 0) break;
        int32_t nf = 0;
        for (int32_t k = 0; k < m; ++k) {
            const int32_t r = odd_tmp_[k];
            for (int32_t x = r; x != -1; x = next_[x]) {
                for (int32_t a = off[x]; a < off[x + 1]; ++a) {
                    const int32_t e = adj[a];
                    if (growth_[e] >= 2) continue;
                    if (find(eu[e]) == find(ev[e])) continue;  // interior edge
                    grow(e, 1, nf);
                }
            }
        }
        for (int32_t k = 0; k < nf; ++k) unite(eu[fused_[k]], ev[fused_[k]]);
        for (int32_t k = 0; k < m; ++k) odd_roots_[k] = odd_tmp_[k];
        n_odd = m;
    }
    // 4. peeling over the spanning forest of fully-grown edges; boundary-rooted tree first
    uint8_t obs = 0;
    for (int32_t pass = 0; pass < 2; ++pass) {
        const int32_t lo = pass == 0 ? B_ : 0, hi = pass == 0 ? B_ + 1 : g_.num_det;
        for (int32_t s = lo; s < hi; ++s) {
            if (visited_[s]) continue;
            if (pass == 1 && !defect_[s]) continue;
            int32_t head = 0, tail = 0;
            queue_[tail++] = s; visited_[s] = 1; mark(s); tree_edge_[s] = -1;
            while (head < tail) {
                const int32_t x = queue_[head++];
                for (int32_t a = off[x]; a < off[x + 1]; ++a) {
                    const int32_t e = adj[a];
                    if (growth_[e] < 2) continue;
                    const int32_t y = eu[e] == x ? ev[e] : eu[e];
                    if (visited_[y]) continue;
                    visited_[y] = 1; mark(y); tree_edge_[y] = e;
                    queue_[tail++] = y;
                }
            }
            for (int32_t k = tail - 1; k > 0; --k) {   // leaves first
                const int32_t y = queue_[k];
                if (!defect_[y]) continue;
                const int32_t e = tree_edge_[y];
                const int32_t x = eu[e] == y ? ev[e] : eu[e];
                defect_[y] = 0; defect_[x] ^= 1; obs ^= g_.eobs[e];
            }
        }
    }
    // 5. reset only what changed (no allocation, O(active))
    for (int32_t k = 0; k < n_grown_; ++k) growth_[grown_list_[k]] = 0;
    for (int32_t k = 0; k < n_dirty_; ++k) {
        const int32_t i = dirty_[k];
        parent_[i] = i; size_[i] = 1; parity_[i] = 0; touches_bnd_[i] = 0;
        defect_[i] = 0; visited_[i] = 0; is_dirty_[i] = 0; next_[i] = -1; tail_[i] = i;
    }
    touches_bnd_[B_] = 1;
    return obs;
}

}  // namespace uf
