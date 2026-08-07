#include "moead.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

using evo::CandidateSolution;
using evo::Matrix;

MOEAD::MOEAD(const Model& model, int pop_size)
    : model_(model), pop_size_(pop_size) {}

std::vector<CandidateSolution> MOEAD::run(const TerminationCriteria& term, ProgressLogger* logger, bool use_repair) {
    std::mt19937 rng(std::random_device{}());

    const auto base_weights = evo::uniform_points_nbi_3d(pop_size_);
    const int n_pop = static_cast<int>(base_weights.size());
    const int T = std::max(2, static_cast<int>(std::ceil(n_pop / 10.0)));

    std::vector<Matrix> prod = evo::random_population(model_, n_pop, rng, 1.5, use_repair);
    std::vector<CandidateSolution> population = evo::evaluate_population(model_, prod);

    std::array<double, 3> zmin = {std::numeric_limits<double>::infinity(),
                                  std::numeric_limits<double>::infinity(),
                                  std::numeric_limits<double>::infinity()};
    for (const auto& c : population) {
        if (c.metrics.unmet_demand > 1e-9) continue;
        zmin[0] = std::min(zmin[0], c.metrics.holding_cost);
        zmin[1] = std::min(zmin[1], c.metrics.setup_cost);
        zmin[2] = std::min(zmin[2], c.metrics.overtime_cost);
    }
    if (!std::isfinite(zmin[0])) {
        for (const auto& c : population) {
            zmin[0] = std::min(zmin[0], c.metrics.holding_cost);
            zmin[1] = std::min(zmin[1], c.metrics.setup_cost);
            zmin[2] = std::min(zmin[2], c.metrics.overtime_cost);
        }
    }

    int generation = 0;
    int evaluations = static_cast<int>(population.size());
    auto start = std::chrono::steady_clock::now();

    while (true) {
        if (term.max_generations > 0 && generation >= term.max_generations) break;
        if (term.max_evaluations > 0 && evaluations >= term.max_evaluations) break;
        if (term.max_seconds > 0.0) {
            double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            if (elapsed >= term.max_seconds) break;
        }
        std::vector<std::vector<int>> neigh(n_pop);
        for (int i = 0; i < n_pop; ++i) {
            std::vector<std::pair<double, int>> dist;
            dist.reserve(n_pop);
            for (int j = 0; j < n_pop; ++j) {
                double d = 0.0;
                for (int k = 0; k < 3; ++k) d += std::pow(base_weights[i][k] - base_weights[j][k], 2);
                dist.push_back({std::sqrt(d), j});
            }
            std::partial_sort(dist.begin(), dist.begin() + std::min(T, static_cast<int>(dist.size())), dist.end(),
                              [](auto& a, auto& b) { return a.first < b.first; });
            for (int k = 0; k < T && k < static_cast<int>(dist.size()); ++k) neigh[i].push_back(dist[k].second);
        }

        for (int i = 0; i < n_pop; ++i) {
            std::vector<int> pool = neigh[i];
            std::shuffle(pool.begin(), pool.end(), rng);

            std::uniform_real_distribution<> ur(0.0, 1.0);
            int a = pool[0];
            int b = pool[1];

            Matrix child = evo::crossover_sbx(prod[a], prod[b], rng, 20.0, evo::ops_moead.pc);
            if (ur(rng) < evo::ops_moead.p_toggle) {
                evo::mutate_toggle_on_off(model_, child, rng);
            }
            if (use_repair && ur(rng) < evo::ops_moead.p_repair) evo::repair_demand(model_, child);
            CandidateSolution cand = evo::evaluate(model_, child);
            evaluations++;

            // Update ideal point only with feasible offspring.
            if (cand.metrics.unmet_demand <= 1e-9) {
                zmin[0] = std::min(zmin[0], cand.metrics.holding_cost);
                zmin[1] = std::min(zmin[1], cand.metrics.setup_cost);
                zmin[2] = std::min(zmin[2], cand.metrics.overtime_cost);
            }

            std::array<double, 3> zmax = {-std::numeric_limits<double>::infinity(),
                                          -std::numeric_limits<double>::infinity(),
                                          -std::numeric_limits<double>::infinity()};
            for (const auto& c : population) {
                if (c.metrics.unmet_demand > 1e-9) continue;
                zmax[0] = std::max(zmax[0], c.metrics.holding_cost);
                zmax[1] = std::max(zmax[1], c.metrics.setup_cost);
                zmax[2] = std::max(zmax[2], c.metrics.overtime_cost);
            }
            if (!std::isfinite(zmax[0])) {
                for (const auto& c : population) {
                    zmax[0] = std::max(zmax[0], c.metrics.holding_cost);
                    zmax[1] = std::max(zmax[1], c.metrics.setup_cost);
                    zmax[2] = std::max(zmax[2], c.metrics.overtime_cost);
                }
            }

            std::array<double, 3> range = {
                std::max(zmax[0] - zmin[0], 1e-9),
                std::max(zmax[1] - zmin[1], 1e-9),
                std::max(zmax[2] - zmin[2], 1e-9)
            };

            // Update neighbours using constrained normalized Tchebycheff.
            for (int idx : pool) {
                auto g_value = [&](const SolutionMetrics& m, const std::array<double, 3>& w) {
                    double v0 = std::abs(m.holding_cost - zmin[0]) / range[0] * w[0];
                    double v1 = std::abs(m.setup_cost - zmin[1]) / range[1] * w[1];
                    double v2 = std::abs(m.overtime_cost - zmin[2]) / range[2] * w[2];
                    return std::max({v0, v1, v2});
                };
                double cv_new = cand.metrics.unmet_demand;
                double cv_old = population[idx].metrics.unmet_demand;
                double g_new = g_value(cand.metrics, base_weights[idx]);
                double g_old = g_value(population[idx].metrics, base_weights[idx]);
                if ((cv_new == cv_old && g_new <= g_old) || (cv_new < cv_old)) {
                    population[idx] = cand;
                    prod[idx] = child;
                }
            }
        }
        generation++;
        if (logger && !population.empty()) {
            const auto* best = &population.front().metrics;
            for (const auto& v : population) {
                if (v.metrics.unmet_demand < best->unmet_demand - 1e-9 ||
                    (std::abs(v.metrics.unmet_demand - best->unmet_demand) < 1e-9 &&
                     v.metrics.total_cost < best->total_cost)) {
                    best = &v.metrics;
                }
            }
            double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            logger->log(generation, evaluations, elapsed, *best);
            logger->log_population_metrics(generation, population);
            logger->log_population(generation, model_, population);
        }
    }

    // Return final population; downstream can filter non-dominated solutions later.
    return population;
}
