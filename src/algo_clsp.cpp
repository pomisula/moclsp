#include "algo_clsp.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <unordered_set>
#include <vector>

#include "evo_ops.h"

namespace {

using Matrix = std::vector<std::vector<double>>;

Solution build_solution(const Model& model, const Matrix& prod) {
    Solution sol(model);
    for (int j = 0; j < model.J; ++j) {
        for (int t = 0; t < model.T; ++t) {
            double qty = prod[j][t];
            if (qty > 0.0) {
                sol.add_production(j, t, 0, qty, 0.0, 0.0);
            }
        }
    }
    return sol;
}

struct MOInd {
    Matrix prod;
    SolutionMetrics metrics;
};

void lns_relax_and_repair(const Model& model, Matrix& prod, std::mt19937& rng, double relax_frac) {
    std::uniform_int_distribution<int> t_dist(0, model.T - 1);
    std::uniform_int_distribution<int> j_dist(0, model.J - 1);
    int relax_periods = std::max(1, std::min(model.T, static_cast<int>(std::round(relax_frac * model.T))));
    std::unordered_set<int> chosen;
    while (static_cast<int>(chosen.size()) < relax_periods) chosen.insert(t_dist(rng));
    for (int t : chosen) {
        for (int j = 0; j < model.J; ++j) prod[j][t] = 0.0;
    }
    for (int k = 0; k < relax_periods; ++k) {
        int j = j_dist(rng);
        int t = t_dist(rng);
        if (t <= 0 || prod[j][t] <= 0.0) continue;
        int s = std::uniform_int_distribution<int>(0, t - 1)(rng);
        double mv = std::max(1.0, prod[j][t] / 2.0);
        prod[j][t] -= mv;
        prod[j][s] += mv;
    }
    evo::repair_demand(model, prod);
}

std::tuple<Matrix, SolutionMetrics, int> lns_search(const Model& model, const Matrix& base, const SolutionMetrics& base_m,
                                                    std::mt19937& rng, int iters, double relax_frac) {
    Matrix best = base;
    SolutionMetrics best_m = base_m;
    int evals = 0;
    for (int i = 0; i < iters; ++i) {
        Matrix tmp = base;
        lns_relax_and_repair(model, tmp, rng, relax_frac);
        SolutionMetrics m = model.evaluate_solution(build_solution(model, tmp));
        ++evals;
        if (evo::dominates(m, best_m)) {
            best = std::move(tmp);
            best_m = m;
        }
    }
    return {std::move(best), best_m, evals};
}

std::vector<int> kmeans_indices(const std::vector<MOInd>& pop, int k, int iters, std::mt19937& rng) {
    int n = static_cast<int>(pop.size());
    if (n == 0 || k <= 0) return {};
    k = std::min(k, n);
    auto normalize = [&](std::vector<std::array<double, 3>>& vals) {
        std::array<double, 3> mn{vals[0][0], vals[0][1], vals[0][2]};
        std::array<double, 3> mx = mn;
        for (const auto& v : vals) {
            for (int d = 0; d < 3; ++d) {
                mn[d] = std::min(mn[d], v[d]);
                mx[d] = std::max(mx[d], v[d]);
            }
        }
        for (auto& v : vals) {
            for (int d = 0; d < 3; ++d) {
                double den = (mx[d] - mn[d]);
                if (den > 0) v[d] = (v[d] - mn[d]) / den;
                else v[d] = 0.0;
            }
        }
    };

    std::vector<std::array<double, 3>> vals(n);
    for (int i = 0; i < n; ++i) {
        const auto& m = pop[i].metrics;
        vals[i] = {m.holding_cost, m.setup_cost, m.overtime_cost};
    }
    normalize(vals);

    std::vector<int> centers;
    centers.reserve(k);
    std::uniform_int_distribution<int> dist_n(0, n - 1);
    centers.push_back(dist_n(rng));
    auto l2 = [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
        double s = 0.0;
        for (int d = 0; d < 3; ++d) {
            double dx = a[d] - b[d];
            s += dx * dx;
        }
        return s;
    };
    while (static_cast<int>(centers.size()) < k) {
        std::vector<double> dist_sq(n, std::numeric_limits<double>::infinity());
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            double best = std::numeric_limits<double>::infinity();
            for (int c : centers) {
                best = std::min(best, l2(vals[i], vals[c]));
            }
            dist_sq[i] = best;
            sum += best;
        }
        if (sum <= 0) break;
        std::uniform_real_distribution<double> ur(0.0, sum);
        double r = ur(rng);
        double acc = 0.0;
        for (int i = 0; i < n; ++i) {
            acc += dist_sq[i];
            if (acc >= r) {
                centers.push_back(i);
                break;
            }
        }
    }

    std::vector<int> assign(n, 0);
    double tol = 1e-6;
    for (int it = 0; it < iters; ++it) {
        for (int i = 0; i < n; ++i) {
            double best = std::numeric_limits<double>::infinity();
            int bestc = 0;
            for (int c = 0; c < k; ++c) {
                double d = l2(vals[i], vals[centers[c]]);
                if (d < best) {
                    best = d;
                    bestc = c;
                }
            }
            assign[i] = bestc;
        }
        std::vector<std::array<double, 3>> sum(k, {0, 0, 0});
        std::vector<int> cnt(k, 0);
        for (int i = 0; i < n; ++i) {
            int c = assign[i];
            for (int d = 0; d < 3; ++d) sum[c][d] += vals[i][d];
            cnt[c]++;
        }
        bool converged = true;
        for (int c = 0; c < k; ++c) {
            if (cnt[c] == 0) continue;
            for (int d = 0; d < 3; ++d) sum[c][d] /= cnt[c];
            double best = std::numeric_limits<double>::infinity();
            int best_idx = centers[c];
            for (int i = 0; i < n; ++i) {
                if (assign[i] != c) continue;
                double d = l2(vals[i], sum[c]);
                if (d < best) {
                    best = d;
                    best_idx = i;
                }
            }
            if (l2(vals[centers[c]], vals[best_idx]) > tol) converged = false;
            centers[c] = best_idx;
        }
        if (converged) break;
    }
    return centers;
}

}  // namespace

AlgoCLSP::AlgoCLSP(const Model& model, int pop_size) : model_(model), pop_size_(pop_size) {}

std::vector<evo::CandidateSolution> AlgoCLSP::run(const TerminationCriteria& term, ProgressLogger* logger) {
    std::mt19937 rng(std::random_device{}());
    std::vector<MOInd> pop;
    pop.reserve(pop_size_);
    int evaluations = 0;

    auto init_prod = evo::random_population(model_, pop_size_, rng, 1.5, true);
    for (auto& p : init_prod) {
        pop.push_back({p, model_.evaluate_solution(build_solution(model_, p))});
        evaluations++;
    }

    int generation = 0;
    auto start = std::chrono::steady_clock::now();
    while (true) {
        if (term.max_generations > 0 && generation >= term.max_generations) break;
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - start).count();
        if (term.max_seconds > 0.0 && elapsed >= term.max_seconds) break;
        if (term.max_evaluations > 0 && evaluations >= term.max_evaluations) break;

        std::vector<MOInd> next = pop;
        while (static_cast<int>(next.size()) < 2 * pop_size_) {
            int a = std::uniform_int_distribution<int>(0, pop_size_ - 1)(rng);
            int b = std::uniform_int_distribution<int>(0, pop_size_ - 1)(rng);
            std::uniform_real_distribution<> ur(0.0, 1.0);
            Matrix child = evo::crossover_sbx(pop[a].prod, pop[b].prod, rng, 20.0, evo::ops_clsp.pc);
            if (ur(rng) < evo::ops_clsp.p_toggle) evo::mutate_toggle_on_off(model_, child, rng);
            // if (ur(rng) < evo::ops_clsp.p_poly_mut) evo::mutate_poly(child, rng, 20.0, evo::ops_clsp.p_poly_mut);
            if (ur(rng) < evo::ops_clsp.p_shift) evo::mutate_shift(model_, child, rng);
            if (ur(rng) < evo::ops_clsp.p_repair) evo::repair_demand(model_, child);
            SolutionMetrics m = model_.evaluate_solution(build_solution(model_, child));
            next.push_back({child, m});
            evaluations += 1;
        }
        auto fronts = evo::fast_non_dominated_sort(next);
        std::vector<MOInd> filtered;
        filtered.reserve(pop_size_);
        for (size_t f = 0; f < fronts.size() && static_cast<int>(filtered.size()) < pop_size_; ++f) {
            auto cd = evo::crowding_distance(next, fronts[f]);
            std::vector<int> order(fronts[f].size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](int a, int b) { return cd[a] > cd[b]; });
            for (int idx_local : order) {
                int idx = fronts[f][idx_local];
                if (static_cast<int>(filtered.size()) >= pop_size_) break;
                filtered.push_back(next[idx]);
            }
        }
        pop.swap(filtered);
        generation++;

        {
            int k = evo::ops_clsp.lns_clusters;
            if (k > 0) {
                auto reps = kmeans_indices(pop, k, 5, rng);
                for (int idx : reps) {
                    auto [lns_prod, lns_metrics, evals] = lns_search(
                        model_, pop[idx].prod, pop[idx].metrics, rng, evo::ops_clsp.lns_iters, evo::ops_clsp.lns_relax_frac);
                    evaluations += evals;
                    if (evo::dominates(lns_metrics, pop[idx].metrics) ) {
                        pop[idx].prod = std::move(lns_prod);
                        pop[idx].metrics = lns_metrics;
                    }
                }
            }
        }

        if (logger) {
            auto it = std::min_element(pop.begin(), pop.end(), [](const MOInd& a, const MOInd& b) {
                return a.metrics.total_cost < b.metrics.total_cost;
            });
            if (it != pop.end()) {
                logger->log(generation, evaluations, elapsed, it->metrics);
                logger->log_population_metrics(generation, pop);
                logger->log_population(generation, model_, pop);
            }
        }
    }

    std::vector<evo::CandidateSolution> result;
    result.reserve(pop.size());
    for (const auto& ind : pop) {
        evo::CandidateSolution res(model_);
        res.plan = build_solution(model_, ind.prod);
        res.metrics = ind.metrics;
        result.push_back(std::move(res));
    }
    return result;
}
