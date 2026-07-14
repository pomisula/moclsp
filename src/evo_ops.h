#pragma once

#include <algorithm>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

#include "model.h"
#include "solution.h"

namespace evo {

using Matrix = std::vector<std::vector<double>>;

struct CandidateSolution {
    Solution plan;
    SolutionMetrics metrics;
    explicit CandidateSolution(const Model& model) : plan(model) {}
};

struct OperatorProbs {
    double pc = 1.0;
    double p_shift = 1.0;
    double p_toggle = 0.0;
    double p_repair = 1.0;
    double p_poly_mut = 0.0;
    int lns_iters = 3;
    double lns_relax_frac = 0.25;
    int lns_clusters = 0;
};

extern OperatorProbs ops_nsgaii;
extern OperatorProbs ops_nsgaiii;
extern OperatorProbs ops_moead;
extern OperatorProbs ops_clsp;

void repair_demand(const Model& model, Matrix& prod);

std::vector<Matrix> random_population(const Model& model, int pop_size, std::mt19937& rng, double scale = 1.5,
                                      bool use_repair = true);

Solution build_solution(const Model& model, const Matrix& prod);

CandidateSolution evaluate(const Model& model, const Matrix& prod);

std::vector<CandidateSolution> evaluate_population(const Model& model, const std::vector<Matrix>& prod_vec);

Matrix crossover_sbx(const Matrix& a, const Matrix& b, std::mt19937& rng, double eta_c = 20.0, double pc = 1.0);

void mutate_shift(const Model& model, Matrix& prod, std::mt19937& rng);
void mutate_poly(Matrix& prod, std::mt19937& rng, double eta_m = 20.0, double pm = -1.0);
void mutate_toggle_on_off(const Model& model, Matrix& prod, std::mt19937& rng,
                          double p_on = -1.0, double p_off = -1.0, double scale = 1.5);

bool dominates(const SolutionMetrics& a, const SolutionMetrics& b, double tol = 1e-9);
std::vector<int> non_dominated_indices(const std::vector<SolutionMetrics>& ms, double tol = 1e-9);

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
