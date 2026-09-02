#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>

using namespace std;

struct Instance {
    int n;
    vector<vector<int>> adj;
    vector<pair<int,int>> edges;
};

static Instance build_instance(long long a, long long b, int N) {
    set<int> distances;
    for (int d = -N; d <= N; ++d) {
        if (d == 0) continue;
        long long v = a * 1LL * d * d + b * 1LL * d;
        v = llabs(v);
        if (0 < v && v < N) distances.insert((int)v);
    }

    Instance I;
    I.n = N;
    I.adj.assign(N, {});
    for (int s : distances) {
        for (int u = 0; u + s < N; ++u) {
            int v = u + s;
            I.edges.push_back({u, v});
            I.adj[u].push_back(v);
            I.adj[v].push_back(u);
        }
    }
    return I;
}

static int energy(const Instance& I, const vector<int>& color) {
    int e = 0;
    for (auto [u,v] : I.edges) e += (color[u] == color[v]);
    return e;
}

static int recolor_delta(const Instance& I, const vector<int>& color,
                         int v, int new_color) {
    int old_color = color[v];
    int delta = 0;
    for (int u : I.adj[v]) {
        if (color[u] == old_color) --delta;
        if (color[u] == new_color) ++delta;
    }
    return delta;
}

static bool valid(const Instance& I, const vector<int>& color) {
    for (auto [u,v] : I.edges)
        if (color[u] == color[v]) return false;
    return true;
}

static void print_word(const vector<int>& color) {
    static const char name[4] = {'R','G','B','Y'};
    for (int c : color) cout << name[c];
    cout << '\n';
}

int main(int argc, char** argv) {
    if (argc < 4) {
        cerr << "usage: " << argv[0]
             << " a b N [seed] [restarts] [steps_per_restart]\n";
        cerr << "example: " << argv[0] << " 1 1 96 1\n";
        return 2;
    }

    long long a = stoll(argv[1]);
    long long b = stoll(argv[2]);
    int N = stoi(argv[3]);
    uint64_t seed = argc >= 5 ? stoull(argv[4]) : random_device{}();
    int restarts = argc >= 6 ? stoi(argv[5]) : 200;
    int steps_per_restart = argc >= 7 ? stoi(argv[6]) : 500000;

    if (N <= 0) {
        cerr << "N must be positive\n";
        return 2;
    }

    Instance I = build_instance(a, b, N);
    mt19937_64 rng(seed);
    uniform_real_distribution<double> U(0.0, 1.0);
    uniform_int_distribution<int> vertex_dist(0, N - 1);
    uniform_int_distribution<int> offset_dist(1, 3);
    uniform_int_distribution<int> color_dist(0, 3);

    vector<int> best_color(N, 0);
    int best_energy = (int)I.edges.size() + 1;

    // The state always has exactly one color per vertex. Therefore the
    // one-hot QUBO penalty is identically zero, and the energy is simply
    // the number of monochromatic forbidden pairs.
    for (int restart = 0; restart < restarts; ++restart) {
        vector<int> color(N);
        for (int& c : color) c = color_dist(rng);
        int E = energy(I, color);

        if (E < best_energy) {
            best_energy = E;
            best_color = color;
        }
        if (E == 0) {
            cout << "found witness: a=" << a << " b=" << b << " N=" << N
                 << " seed=" << seed << " restart=" << restart << "\n";
            print_word(color);
            return valid(I, color) ? 0 : 1;
        }

        double T0 = max(0.5, 0.25 * sqrt(max(1.0, 2.0 * I.edges.size() / N)));
        double Tmin = 1e-3;
        double alpha = pow(Tmin / T0, 1.0 / max(1, steps_per_restart));
        double T = T0;

        for (int step = 0; step < steps_per_restart; ++step) {
            int v = vertex_dist(rng);
            int old = color[v];
            int q = (old + offset_dist(rng)) & 3;
            int delta = recolor_delta(I, color, v, q);

            if (delta <= 0 || U(rng) < exp(-double(delta) / T)) {
                color[v] = q;
                E += delta;

                if (E < best_energy) {
                    best_energy = E;
                    best_color = color;
                }
                if (E == 0) {
                    cout << "found witness: a=" << a << " b=" << b
                         << " N=" << N << " seed=" << seed
                         << " restart=" << restart << " step=" << step << "\n";
                    print_word(color);
                    return valid(I, color) ? 0 : 1;
                }
            }
            T *= alpha;
        }
    }

    cerr << "no zero-energy coloring found; best energy=" << best_energy
         << ", seed=" << seed << "\n";
    print_word(best_color);
    return 1;
}
