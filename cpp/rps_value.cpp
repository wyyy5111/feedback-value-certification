// SPDX-License-Identifier: MIT
#include <algorithm>
#include <array>
#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Deterministic interval arithmetic for the declared ideal feedback model.
// The state and finite feedback kernels use exact rational parameters.
//
// y_i=(1/50)c_i+(1-s_t)/3, sum(y)=1, y_0=(1/3,1/3,1/3).
// y'=(49/50)y+(1/50)e_b. Adding the common shift changes neither selector.
// Coordinates below are u=G*y_R, v=G*y_P. A piece is a CLOSED grid box
// intersected with the simplex and one weak principal-action region. Its tag
// fixes the current principal action. Ties are an outer relaxation, not a
// claim that every tag equals the source's deterministic tie-breaking rule.
// The exact initial tie is evaluated separately with its source action.
//
// An edge exists only when the source piece intersects the inverse image of
// the target piece under F_b. The SAME target piece is used for all next-WL
// states j; no independent minimization over the three image cells occurs.

namespace {
using I = __int128_t;
using U = __uint128_t;
constexpr std::int64_t Q = std::int64_t{1} << 40;
constexpr std::uint64_t OWNED_BYTE_GUARD = 1024ULL * 1024 * 1024;

struct Guard : std::runtime_error { using std::runtime_error::runtime_error; };
struct HP { std::int64_t a, b, c; }; // a*u+b*v <= c
struct Affine { std::int64_t a, b, c; };
struct TargetCells { std::array<std::uint32_t, 4> ids{}; std::uint8_t count = 0; };
struct Images { std::array<TargetCells, 3> by_b; };
struct Piece { std::uint32_t cell; std::uint8_t principal; };

struct GeometryCounts {
    std::uint64_t feasibility_calls = 0;
    std::uint64_t corner_tests = 0;
    std::uint64_t intersection_pairs = 0;
};

struct Graph {
    int G;
    std::string selector;
    std::vector<Piece> pieces;
    std::vector<std::array<std::int32_t, 3>> at_cell;
    std::vector<std::uint32_t> offsets, edges;
    std::array<std::array<HP, 2>, 3> cuts;
    GeometryCounts geometry;
};

struct Value {
    std::int64_t lo = 0, hi = 0;
    std::uint64_t pieces = 0, edges = 0, processed_states = 0;
    std::uint64_t candidate_rows = 0, kernel_terms = 0, weighted_products = 0, outward_divisions = 0;
    GeometryCounts geometry;
};

std::string integer_string(I n) {
    if (n == 0) return "0";
    const bool negative = n < 0;
    U x = negative ? U(-(n + 1)) + 1 : U(n);
    std::string s;
    while (x) { s.push_back(char('0' + unsigned(x % 10))); x /= 10; }
    if (negative) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}

std::int64_t narrow(I n) {
    if (n < std::numeric_limits<std::int64_t>::min()
            || n > std::numeric_limits<std::int64_t>::max())
        throw Guard("Fixed-point value exceeds int64 storage");
    return std::int64_t(n);
}

std::int64_t floor_div(I n, std::int64_t d) {
    I z = n / d;
    if (n % d < 0) --z;
    return narrow(z);
}

std::int64_t ceil_div(I n, std::int64_t d) {
    I z = n / d;
    if (n % d > 0) ++z;
    return narrow(z);
}

std::uint64_t counter(U n) {
    if (n > std::numeric_limits<std::uint64_t>::max())
        throw Guard("Analytical operation counter would overflow");
    return std::uint64_t(n);
}

int score(int a, int b) {
    const int d = (b - a + 3) % 3;
    return d == 2 ? 1 : d == 1 ? -1 : 0;
}

std::array<int, 3> kernel(int principal, int b) {
    std::array<int, 3> s = {1, 1, 1}, k{};
    s[principal] = 58;
    for (int j = 0; j < 3; ++j) {
        const int a = (j + 2) % 3;
        k[j] = s[b] + (a != b ? 3 * s[a] : 0);
    }
    if (k[0] + k[1] + k[2] != 180)
        throw std::logic_error("WL kernel does not sum to 180");
    return k;
}

std::array<HP, 2> principal_cuts(const std::string& selector, int h, int G) {
    std::array<Affine, 3> f;
    int chosen;
    if (selector == "mode") {
        f = {{{1, 0, 0}, {0, 1, 0}, {-1, -1, G}}};
        chosen = (h + 2) % 3; // h beats the modal predicted opponent action.
    } else if (selector == "BR") {
        f = {{{-1, -2, G}, {2, 1, -G}, {-1, 1, 0}}};
        chosen = h;
    } else throw std::logic_error("Unknown selector");
    std::array<HP, 2> result{};
    int n = 0;
    for (int j = 0; j < 3; ++j) if (j != chosen)
        result[n++] = {f[j].a - f[chosen].a, f[j].b - f[chosen].b,
                       f[chosen].c - f[j].c};
    return result;
}

// Axis box coordinates are rational with the common POSITIVE denominator den.
// The remaining <=5 halfplanes are exact integers. Fast corner checks usually
// settle interior cells. Otherwise a bounded closed nonempty polygon has a
// vertex given by two nonparallel active boundary lines, including degeneracies.
bool feasible_box(std::int64_t lu, std::int64_t hu,
                  std::int64_t lv, std::int64_t hv, std::int64_t den,
                  const std::array<HP, 5>& cuts, int nc, GeometryCounts& counts) {
    ++counts.feasibility_calls;
    if (lu > hu || lv > hv) return false;
    for (std::int64_t u : {lu, hu}) for (std::int64_t v : {lv, hv}) {
        ++counts.corner_tests;
        bool good = true;
        for (int k = 0; k < nc; ++k)
            if (I(cuts[k].a) * u + I(cuts[k].b) * v > I(cuts[k].c) * den) {
                good = false; break;
            }
        if (good) return true;
    }
    std::array<HP, 9> all{};
    all[0] = {-den, 0, -lu}; all[1] = {den, 0, hu};
    all[2] = {0, -den, -lv}; all[3] = {0, den, hv};
    for (int k = 0; k < nc; ++k) all[4 + k] = cuts[k];
    const int n = 4 + nc;
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) {
        ++counts.intersection_pairs;
        I d = I(all[i].a) * all[j].b - I(all[j].a) * all[i].b;
        if (d == 0) continue;
        I un = I(all[i].c) * all[j].b - I(all[j].c) * all[i].b;
        I vn = I(all[i].a) * all[j].c - I(all[j].a) * all[i].c;
        if (d < 0) { d = -d; un = -un; vn = -vn; }
        bool good = true;
        for (int k = 0; k < n; ++k)
            if (I(all[k].a) * un + I(all[k].b) * vn > I(all[k].c) * d) {
                good = false; break;
            }
        if (good) return true;
    }
    return false;
}

std::pair<int, int> closed_index_range(std::int64_t low, std::int64_t high,
                                      std::int64_t denominator, int G) {
    // All numerators here are nonnegative. If the lower endpoint is exactly
    // on a grid boundary, the preceding closed cell must also be included.
    int first = int(low / denominator);
    if (low % denominator == 0) --first;
    int last = int(high / denominator);
    first = std::max(0, first); last = std::min(G - 1, last);
    return {first, last};
}

void preflight(int T, int G) {
    const U cells = U(G) * G, p = 3 * cells;
    const U e = 36 * p; // 3 b * <=4 image cells * <=3 target pieces.
    const U bytes = cells * (sizeof(Images) + sizeof(std::array<std::int32_t, 3>))
        + p * sizeof(Piece) + (3 * p + 1) * sizeof(std::uint32_t)
        + e * sizeof(std::uint32_t) + 12 * p * sizeof(std::int64_t);
    if (bytes > OWNED_BYTE_GUARD)
        throw Guard("Conservative owned-memory preflight exceeds 1024 MiB");
    if (p > std::numeric_limits<std::int32_t>::max()
            || e > std::numeric_limits<std::uint32_t>::max())
        throw Guard("Piece or edge index exceeds stored integer type");
    if (I(T) * (Q + 1) > std::numeric_limits<std::int64_t>::max() / 2)
        throw Guard("Horizon exceeds fixed-point storage preflight");
    counter(U(T) * e * 6 + 100);
}

std::vector<Images> make_images(int G) {
    const std::size_t cells = std::size_t(G) * G;
    std::vector<Images> images(cells);
    for (int i = 0; i < G; ++i) for (int v = 0; v < G; ++v) {
        const std::size_t cell = std::size_t(i) * G + v;
        for (int b = 0; b < 3; ++b) {
            const std::int64_t ul = 49LL * i + (b == 0 ? G : 0);
            const std::int64_t vl = 49LL * v + (b == 1 ? G : 0);
            const auto ur = closed_index_range(ul, ul + 49, 50, G);
            const auto vr = closed_index_range(vl, vl + 49, 50, G);
            auto& out = images[cell].by_b[b];
            for (int u = ur.first; u <= ur.second; ++u)
                for (int w = vr.first; w <= vr.second; ++w) {
                    if (out.count == out.ids.size())
                        throw std::logic_error("Shrinking box intersects more than four closed cells");
                    out.ids[out.count++] = std::uint32_t(std::size_t(u) * G + w);
                }
            if (!out.count) throw std::logic_error("Empty affine image cover");
        }
    }
    return images;
}

Graph make_graph(int G, const std::string& selector, const std::vector<Images>& images) {
    Graph graph;
    graph.G = G; graph.selector = selector;
    const std::size_t cells = std::size_t(G) * G;
    graph.at_cell.assign(cells, {{-1, -1, -1}});
    graph.pieces.reserve(3 * cells);
    for (int h = 0; h < 3; ++h) graph.cuts[h] = principal_cuts(selector, h, G);
    for (int i = 0; i < G; ++i) for (int j = 0; j < G; ++j) {
        if (i + j > G) continue;
        const std::size_t cell = std::size_t(i) * G + j;
        for (int h = 0; h < 3; ++h) {
            std::array<HP, 5> cuts{};
            cuts[0] = {1, 1, G};
            cuts[1] = graph.cuts[h][0]; cuts[2] = graph.cuts[h][1];
            if (feasible_box(i, i + 1, j, j + 1, 1, cuts, 3, graph.geometry)) {
                graph.at_cell[cell][h] = std::int32_t(graph.pieces.size());
                graph.pieces.push_back({std::uint32_t(cell), std::uint8_t(h)});
            }
        }
    }
    graph.offsets.reserve(3 * graph.pieces.size() + 1);
    graph.edges.reserve(36 * graph.pieces.size());
    graph.offsets.push_back(0);
    for (const Piece& piece : graph.pieces) {
        const int i = int(piece.cell / G), v = int(piece.cell % G);
        for (int b = 0; b < 3; ++b) {
            const std::size_t begin = graph.edges.size();
            const int shift_u = b == 0 ? G : 0, shift_v = b == 1 ? G : 0;
            const auto& targets = images[piece.cell].by_b[b];
            for (int ti = 0; ti < targets.count; ++ti) {
                const std::uint32_t target_cell = targets.ids[ti];
                const int u = int(target_cell / G), w = int(target_cell % G);
                // Combined source/target axis box, denominator 49.
                const std::int64_t lu = std::max(49LL * i, 50LL * u - shift_u);
                const std::int64_t hu = std::min(49LL * (i + 1), 50LL * (u + 1) - shift_u);
                const std::int64_t lv = std::max(49LL * v, 50LL * w - shift_v);
                const std::int64_t hv = std::min(49LL * (v + 1), 50LL * (w + 1) - shift_v);
                if (lu > hu || lv > hv) continue;
                for (int h = 0; h < 3; ++h) {
                    const std::int32_t target_piece = graph.at_cell[target_cell][h];
                    if (target_piece < 0) continue;
                    std::array<HP, 5> cuts{};
                    cuts[0] = {1, 1, G};
                    cuts[1] = graph.cuts[piece.principal][0];
                    cuts[2] = graph.cuts[piece.principal][1];
                    for (int k = 0; k < 2; ++k) {
                        const HP q = graph.cuts[h][k];
                        cuts[3 + k] = {49 * q.a, 49 * q.b,
                            50 * q.c - q.a * shift_u - q.b * shift_v};
                    }
                    // The target simplex inequality is redundant: F_b maps
                    // the source simplex into itself exactly.
                    if (feasible_box(lu, hu, lv, hv, 49, cuts, 5, graph.geometry))
                        graph.edges.push_back(std::uint32_t(target_piece));
                }
            }
            if (graph.edges.size() == begin)
                throw std::logic_error("A nonempty source piece has no covered image piece");
            graph.offsets.push_back(std::uint32_t(graph.edges.size()));
        }
    }
    if (graph.offsets.size() != graph.pieces.size() * 3 + 1)
        throw std::logic_error("CSR edge table size mismatch");
    return graph;
}

bool point_in_principal(const Graph& graph, int h, I un, I vn, I den) {
    for (const HP& q : graph.cuts[h])
        if (I(q.a) * un + I(q.b) * vn > I(q.c) * den) return false;
    return true;
}

int exact_principal(const Graph& graph,I un,I vn,I den) {
    const I wn=I(graph.G)*den-un-vn;
    std::array<I,3> f=graph.selector=="mode" ? std::array<I,3>{un,vn,wn}
                                                     : std::array<I,3>{wn-vn,un-wn,vn-un};
    int selected=0;
    for(int j=1;j<3;++j)if(f[j]>f[selected])selected=j;
    return graph.selector=="mode" ? (selected+1)%3 : selected;
}
std::pair<std::int64_t,std::int64_t> exact_prefix(
        const Graph& graph,const std::vector<std::int64_t>& lo,
        const std::vector<std::int64_t>& hi,I un,I vn,I den,int b,int depth,
        std::uint64_t& prefix_nodes) {
    ++prefix_nodes;
    if(depth==0){
        const int h=exact_principal(graph,un,vn,den);
        const auto ur=closed_index_range(narrow(un),narrow(un),narrow(den),graph.G);
        const auto vr=closed_index_range(narrow(vn),narrow(vn),narrow(den),graph.G);
        std::int64_t lower=std::numeric_limits<std::int64_t>::max();
        std::int64_t upper=std::numeric_limits<std::int64_t>::min();bool any=false;
        for(int u=ur.first;u<=ur.second;++u)for(int v=vr.first;v<=vr.second;++v){
            const auto id=graph.at_cell[std::size_t(u)*graph.G+v][h];
            if(id>=0 && point_in_principal(graph,h,un,vn,den)){
                lower=std::min(lower,lo[std::size_t(id)*3+b]);
                upper=std::max(upper,hi[std::size_t(id)*3+b]);any=true;
            }
        }
        if(!any)throw std::logic_error("Exact prefix leaf has no actual-tag cover");
        return {lower,upper};
    }
    const int h=exact_principal(graph,un,vn,den);const auto k=kernel(h,b);
    const I nu=49*un+(b==0?I(graph.G)*den:0);
    const I nv=49*vn+(b==1?I(graph.G)*den:0);
    I lower=I(171)*Q*score(h,b),upper=lower;
    for(int j=0;j<3;++j){
        const auto c=exact_prefix(graph,lo,hi,nu,nv,50*den,j,depth-1,prefix_nodes);
        lower+=I(k[j])*c.first;upper+=I(k[j])*c.second;
    }
    return {floor_div(lower,180),ceil_div(upper,180)};
}
Value evaluate(int T, const Graph& graph) {
    const std::size_t states = 3 * graph.pieces.size();
    std::vector<std::int64_t> next_lo(states, 0), next_hi(states, 0), lo(states), hi(states);
    std::array<std::array<std::array<int, 3>, 3>, 3> k{};
    for (int h = 0; h < 3; ++h) for (int b = 0; b < 3; ++b) k[h][b] = kernel(h, b);
    const I infinity = I(1) << 120;
    for (int remaining = 1; remaining <= T - std::min(T,8); ++remaining) {
        for (std::size_t p = 0; p < graph.pieces.size(); ++p) {
            const int h = graph.pieces[p].principal;
            for (int b = 0; b < 3; ++b) {
                const std::size_t state = 3 * p + b;
                I minimum = infinity, maximum = -infinity;
                const I reward = I(171) * Q * score(h, b);
                const auto weights = k[h][b];
                for (std::uint32_t e = graph.offsets[state]; e < graph.offsets[state + 1]; ++e) {
                    const std::size_t target = std::size_t(graph.edges[e]) * 3;
                    const I nl = reward + I(weights[0]) * next_lo[target]
                        + I(weights[1]) * next_lo[target + 1] + I(weights[2]) * next_lo[target + 2];
                    const I nh = reward + I(weights[0]) * next_hi[target]
                        + I(weights[1]) * next_hi[target + 1] + I(weights[2]) * next_hi[target + 2];
                    minimum = std::min(minimum, nl); maximum = std::max(maximum, nh);
                }
                // All candidate denominators are 180. Min commutes with floor,
                // and max with ceil, so rounding once after extrema is exactly
                // equivalent to directed rounding of every Bellman row.
                lo[state] = floor_div(minimum, 180); hi[state] = ceil_div(maximum, 180);
                if (lo[state] > hi[state]) throw std::logic_error("Reversed Bellman interval");
            }
        }
        next_lo.swap(lo); next_hi.swap(hi);
        if (remaining % 100 == 0 || remaining == T - std::min(T,8))
            std::cerr << graph.selector << " completed_layers=" << remaining
                      << '/' << T - std::min(T,8) << " pieces=" << graph.pieces.size()
                      << " edges=" << graph.edges.size() << '\n';
    }
    Value result;

    I sum_lo=0,sum_hi=0;std::uint64_t prefix_nodes=0;
    for(int b=0;b<3;++b){
        const auto bound=exact_prefix(graph,next_lo,next_hi,graph.G,graph.G,3,b,std::min(T,8),prefix_nodes);
        sum_lo+=bound.first;sum_hi+=bound.second;
    }
    result.lo=floor_div(sum_lo,3);result.hi=ceil_div(sum_hi,3);
    const std::uint64_t initial_rows=prefix_nodes;
    result.pieces = graph.pieces.size(); result.edges = graph.edges.size();
    result.processed_states = counter(U(T - std::min(T,8)) * states + initial_rows);
    result.candidate_rows = counter(U(T - std::min(T,8)) * graph.edges.size() + initial_rows);
    result.kernel_terms = counter(U(3) * result.candidate_rows);
    result.weighted_products = counter(U(6) * result.candidate_rows);
    result.outward_divisions = counter(U(2) * result.processed_states + 2);
    result.geometry = graph.geometry;
    return result;
}

void print_value(const char* name, const Value& v) {
    std::cout << '"' << name << "\":{\"lower_numerator\":" << v.lo
              << ",\"upper_numerator\":" << v.hi << ",\"denominator\":" << Q
              << ",\"pieces\":" << v.pieces << ",\"edges\":" << v.edges
              << ",\"processed_bellman_states\":" << v.processed_states
              << ",\"candidate_rows\":" << v.candidate_rows
              << ",\"analytical_kernel_terms\":" << v.kernel_terms
              << ",\"weighted_fixedpoint_products\":" << v.weighted_products
              << ",\"outward_integer_divisions\":" << v.outward_divisions
              << ",\"geometry_feasibility_calls\":" << v.geometry.feasibility_calls
              << ",\"geometry_corner_tests\":" << v.geometry.corner_tests
              << ",\"geometry_vertex_pairs\":" << v.geometry.intersection_pairs << '}';
}

int positive_int(const std::string& text) {
    std::size_t used = 0;
    const unsigned long long x = std::stoull(text, &used);
    if (used != text.size() || text.empty() || text.front() == '-'
            || x == 0 || x > unsigned(std::numeric_limits<int>::max()))
        throw std::invalid_argument("Positive integer parameter required");
    return int(x);
}
} // namespace

int main(int argc, char** argv) {
    try {
        const char* usage = "Usage: rps-value --T HORIZON --G GRID_RESOLUTION\n"
                            "Both parameters must be positive integers. The exact prefix is min(T, 8).\n";
        if (argc == 1) { std::cout << usage; return 0; }
        int T = 0, G = 0;
        for (int i = 1; i < argc; ++i) {
            const std::string option = argv[i];
            if (option == "--help" || option == "-h") { std::cout << usage; return 0; }
            if ((option == "--T" || option == "--G") && i + 1 < argc) {
                const int x = positive_int(argv[++i]);
                if (option == "--T") T = x; else G = x;
            } else throw std::invalid_argument(usage);
        }
        if (T == 0 || G == 0) throw std::invalid_argument(usage);
        preflight(T, G);
        const auto images = make_images(G);
        Value mode, br;
        { const Graph graph = make_graph(G, "mode", images); mode = evaluate(T, graph); }
        { const Graph graph = make_graph(G, "BR", images); br = evaluate(T, graph); }
        std::cout << "{\"status\":\"COMPLETE_EXACT_PREFIX_STRENGTHENED_GRID_INTERVAL\","
                  << "\"T\":" << T << ",\"G\":" << G << ",\"grid_width_denominator\":" << G
                  << ",\"fixed_point_Q\":" << Q << ",\"exact_prefix_rounds\":" << std::min(T,8) << ",\"grid_layers\":" << T - std::min(T,8) << ',';
        print_value("mode", mode); std::cout << ','; print_value("BR", br);
        std::cout << ",\"BR_minus_mode\":{\"lower_numerator\":"
                  << integer_string(I(br.lo) - mode.hi) << ",\"upper_numerator\":"
                  << integer_string(I(br.hi) - mode.lo) << ",\"denominator\":" << Q << '}'
                  << ",\"total_processed_bellman_states\":"
                  << counter(U(mode.processed_states) + br.processed_states)
                  << ",\"total_analytical_kernel_terms\":"
                  << counter(U(mode.kernel_terms) + br.kernel_terms)
                  << ",\"total_weighted_fixedpoint_products\":"
                  << counter(U(mode.weighted_products) + br.weighted_products)
                  << ",\"initial_principal\":{\"mode\":1,\"BR\":0},"
                  << "\"initial_state\":\"exact_uniform_y_and_uniform_WL_b\","
                  << "\"geometry\":\"closed_cell_simplex_weak_principal_pieces_exact_integer_feasibility\","
                  << "\"same_target_piece_for_all_next_WL_states\":true,"
                  << "\"owned_memory_preflight_guard_bytes\":" << OWNED_BYTE_GUARD << ','
                  << "\"arithmetic\":\"directed_integer_rounding\"}\n";
        return 0;
    } catch (const Guard& e) {
        std::cerr << e.what() << '\n';
        std::cout << "{\"status\":\"UNKNOWN_RESOURCE_PREFLIGHT_GUARD\"}\n";
        return 2;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        std::cout << "{\"status\":\"UNKNOWN_IMPLEMENTATION_OR_ALLOCATION_FAILURE\"}\n";
        return 1;
    }
}
