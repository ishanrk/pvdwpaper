#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
#include <vector>
using namespace std;

struct Solver {
    int n;
    vector<vector<int>> adj;
    vector<int> color, uncolored_degree;
    vector<array<unsigned short,4>> neighbor_color_count;
    vector<unsigned char> neighbor_color_mask;
    uint64_t nodes = 0;

    explicit Solver(vector<vector<int>> graph)
      : n((int)graph.size()), adj(move(graph)), color(n,-1),
        uncolored_degree(n), neighbor_color_count(n), neighbor_color_mask(n,0) {
        for (int v=0; v<n; ++v) {
            uncolored_degree[v] = (int)adj[v].size();
            neighbor_color_count[v].fill(0);
        }
    }

    int choose_vertex() const {
        int best=-1, best_sat=-1, best_deg=-1;
        for (int v=0; v<n; ++v) if (color[v] < 0) {
            int sat = __builtin_popcount((unsigned)neighbor_color_mask[v]);
            if (sat > best_sat || (sat == best_sat && uncolored_degree[v] > best_deg)) {
                best=v; best_sat=sat; best_deg=uncolored_degree[v];
            }
        }
        return best;
    }

    bool assign(int v, int c) {
        color[v]=c;
        bool ok=true;
        for (int u: adj[v]) if (color[u] < 0) {
            --uncolored_degree[u];
            if (neighbor_color_count[u][c]++ == 0)
                neighbor_color_mask[u] |= (1u<<c);
            if (neighbor_color_mask[u] == 15) ok=false;
        }
        return ok;
    }

    void unassign(int v, int c) {
        for (int u: adj[v]) if (color[u] < 0) {
            if (--neighbor_color_count[u][c] == 0)
                neighbor_color_mask[u] &= ~(1u<<c);
            ++uncolored_degree[u];
        }
        color[v]=-1;
    }

    bool dfs(int max_used, int left) {
        ++nodes;
        if (left == 0) return true;
        int v = choose_vertex();
        unsigned forbidden = neighbor_color_mask[v];
        int top = min(3, max_used+1); // existing colors plus at most one new color
        for (int c=0; c<=top; ++c) if (!(forbidden & (1u<<c))) {
            bool is_new = (c == max_used+1);
            bool ok = assign(v,c);
            if (ok && dfs(is_new ? c : max_used, left-1)) return true;
            unassign(v,c);
        }
        return false;
    }
};

static void add_edge(vector<vector<int>>& g, int a, int b) {
    g[a].push_back(b);
    g[b].push_back(a);
}

static int grid_id(int x, int y, int Y) { return x*(Y+1)+y; }

vector<vector<int>> build_grid_graph_even(int X=10, int Y=64) {
    vector<pair<int,int>> pts;
    vector<vector<int>> id(X+1, vector<int>(Y+1,-1));
    for (int x=0; x<=X; ++x)
        for (int y=0; y<=Y; ++y)
            if (((x+y)&1)==0) {
                id[x][y]=(int)pts.size();
                pts.push_back({x,y});
            }
    vector<vector<int>> g(pts.size());
    for (int d=1; d<=X && d*d<=Y; ++d) {
        for (int x=0; x+d<=X; ++x) {
            for (int y=0; y<=Y; ++y) if (id[x][y]>=0) {
                if (y+d*d<=Y && id[x+d][y+d*d]>=0)
                    add_edge(g,id[x][y],id[x+d][y+d*d]);
                if (y-d*d>=0 && id[x+d][y-d*d]>=0)
                    add_edge(g,id[x][y],id[x+d][y-d*d]);
            }
        }
    }
    return g;
}

vector<vector<int>> build_polynomial_graph(int a, int b, int N) {
    set<int> distances;
    // For the small instances checked here, |d| <= N is far more than sufficient.
    for (int d=-N; d<=N; ++d) if (d != 0) {
        long long v = 1LL*a*d*d + 1LL*b*d;
        if (v < 0) v = -v;
        if (v > 0 && v < N) distances.insert((int)v);
    }
    vector<vector<int>> g(N);
    for (int s: distances)
        for (int i=0; i+s<N; ++i)
            add_edge(g, i, i+s);
    return g;
}

bool verify_unsat(vector<vector<int>> g, const string& name) {
    long long twice_edges=0;
    for (const auto& a: g) twice_edges += (long long)a.size();
    Solver s(move(g));
    bool sat = s.dfs(-1, s.n);
    cout << name << ": vertices=" << s.n << " edges=" << twice_edges/2
         << " result=" << (sat ? "SAT" : "UNSAT") << "\n";
    return !sat;
}

int main(int argc, char** argv) {
    string mode = argc >= 2 ? argv[1] : "small";
    bool ok=true;

    if (mode == "small" || mode == "all") {
        for (int r=2; r<=8; ++r)
            ok &= verify_unsat(build_polynomial_graph(1,r,97),
                               "W(x^2+" + to_string(r) + "x;4) <= 97");
    }

    if (mode == "grid" || mode == "all")
        ok &= verify_unsat(build_grid_graph_even(), "Q_even[0..10]x[0..64]");

    if (!ok) {
        cerr << "At least one required UNSAT claim failed.\n";
        return 1;
    }
    return 0;
}
