#pragma once

#include <algorithm>
#include <array>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

#include "model.h"
#include "solution.h"

// Common evolutionary helper types/operators so that algorithms stay "clean".
// Swap/disable operators here when experimenting.
namespace evo {

using Matrix = std::vector<std::vector<double>>;

struct CandidateSolution {
    Solution plan;
    SolutionMetrics metrics;
    explicit CandidateSolution(const Model& model) : plan(model) {}
};

struct OperatorProbs {
    double pc = 1.0;         // SBX crossover probability (per cell inside SBX)
    double p_shift = 1.0;    // outer shift trigger
    double p_toggle = 0.0;   // outer toggle trigger; per-cell prob inside operator
    double p_repair = 1.0;
    double p_lns = 0.0;
    int lns_iters = 3;
    double lns_relax_frac = 0.25;  // fraction of periods to relax in LNS
    int lns_clusters = 0;          // number of clusters for LNS refinement (0 => pop_size/5)
};

// Per-algorithm operator probability settings (configured in main.cpp).
extern OperatorProbs ops_nsgaii;
extern OperatorProbs ops_nsgaiii;
extern OperatorProbs ops_moead;
extern OperatorProbs ops_clsp;

// Convert demands to a dense matrix (double).
Matrix demand_matrix(const Model& model);

// Ensure all demand is covered (adds production if needed), continuous.
void repair_demand(const Model& model, Matrix& prod);

// Random production matrices U(0, scale * demand); optional repair.
std::vector<Matrix> random_population(const Model& model, int pop_size, std::mt19937& rng, double scale = 1.5,
                                      bool use_repair = true);

// Build a Solution object from a production plan.
Solution build_solution(const Model& model, const Matrix& prod);

// Evaluate a production plan and return metrics.
CandidateSolution evaluate(const Model& model, const Matrix& prod);

// Evaluate many production plans.
std::vector<CandidateSolution> evaluate_population(const Model& model, const std::vector<Matrix>& prod_vec);

// PlatEMO-style UniformPoint(N,3,'NBI'); returned size is the actual sample count.
std::vector<std::array<double, 3>> uniform_points_nbi_3d(int n);

// Basic genetic operators (problem-aware). You may comment/replace.
Matrix crossover_mix(const Matrix& a, const Matrix& b, std::mt19937& rng);
// Simulated binary crossover (SBX) on all cells.
Matrix crossover_sbx(const Matrix& a, const Matrix& b, std::mt19937& rng, double eta_c = 20.0, double pc = 1.0);
// Cell-level uniform crossover (per (j,t) choose parent with prob p_b)
Matrix crossover_uniform_cells(const Matrix& a, const Matrix& b, std::mt19937& rng, double p_b = 0.5);

void mutate_shift(const Model& model, Matrix& prod, std::mt19937& rng);
// Polynomial mutation on all cells (continuous).
void mutate_poly(Matrix& prod, std::mt19937& rng, double eta_m = 20.0, double pm = -1.0);
void lns_relax(const Model& model, Matrix& prod, std::mt19937& rng, int relax_periods = -1);

// Sparsify: randomly set a fraction of cells to zero (helps sparse 0/1 style)
void mutate_zero_sparse(const Model& model, Matrix& prod, std::mt19937& rng, double drop_prob = -1.0);

// Local resample: with probability p reset some cells to U(0, scale*demand)
void mutate_resample(const Model& model, Matrix& prod, std::mt19937& rng, double resample_prob = -1.0,
                     double scale = 1.5);

// Explicit on/off flip: if zero turn on with p_on; if >0 switch off with p_off
void mutate_toggle_on_off(const Model& model, Matrix& prod, std::mt19937& rng,
                          double p_on = -1.0, double p_off = -1.0, double scale = 1.5);

// Tournament selection over a score vector (lower is better).
int tournament_select(const std::vector<double>& score, int k, std::mt19937& rng);

// Pareto dominance utilities
bool dominates(const SolutionMetrics& a, const SolutionMetrics& b, double tol = 1e-9);
std::vector<int> non_dominated_indices(const std::vector<SolutionMetrics>& ms, double tol = 1e-9);

// Crowding distance for one front; pop[i].metrics is accessed.
template <typename Vec>
std::vector<double> crowding_distance(const Vec& pop, const std::vector<int>& front) {
    int m = static_cast<int>(front.size());
    std::vector<double> dist(m, 0.0);
    if (m == 0) return dist;
    auto add_obj = [&](auto getter) {
        std::vector<int> idx(m);
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b) {
            return getter(pop[front[a]].metrics) < getter(pop[front[b]].metrics);
        });
        double minv = getter(pop[front[idx.front()]].metrics);
        double maxv = getter(pop[front[idx.back()]].metrics);
        dist[idx.front()] = dist[idx.back()] = std::numeric_limits<double>::infinity();
        if (maxv - minv <= 0) return;
        for (int k = 1; k < m - 1; ++k) {
            double prev = getter(pop[front[idx[k - 1]]].metrics);
            double next = getter(pop[front[idx[k + 1]]].metrics);
            dist[idx[k]] += (next - prev) / (maxv - minv);
        }
    };
    add_obj([](const SolutionMetrics& s) { return s.holding_cost; });
    add_obj([](const SolutionMetrics& s) { return s.setup_cost; });
    add_obj([](const SolutionMetrics& s) { return s.overtime_cost; });
    add_obj([](const SolutionMetrics& s) { return s.unmet_demand; });
    return dist;
}

// Fast non-dominated sort (fronts of indices) for a population exposing .metrics
template <typename Vec>
std::vector<std::vector<int>> fast_non_dominated_sort(const Vec& pop, double tol = 1e-9) {
    int n = static_cast<int>(pop.size());
    std::vector<int> dom_count(n, 0);
    std::vector<std::vector<int>> S(n);
    std::vector<std::vector<int>> fronts;
    for (int p = 0; p < n; ++p) {
        for (int q = 0; q < n; ++q) {
            if (p == q) continue;
            if (dominates(pop[p].metrics, pop[q].metrics, tol)) {
                S[p].push_back(q);
            } else if (dominates(pop[q].metrics, pop[p].metrics, tol)) {
                dom_count[p]++;
            }
        }
        if (dom_count[p] == 0) {
            if (fronts.empty()) fronts.push_back({});
            fronts[0].push_back(p);
        }
    }
    int i = 0;
    while (i < static_cast<int>(fronts.size()) && !fronts[i].empty()) {
        std::vector<int> next;
        for (int p : fronts[i]) {
            for (int q : S[p]) {
                dom_count[q]--;
                if (dom_count[q] == 0) next.push_back(q);
            }
        }
        if (!next.empty()) fronts.push_back(next);
        ++i;
    }
    return fronts;
}

}  // namespace evo
