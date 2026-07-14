#include "nsgaii.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>
#include <chrono>
#include <numeric>

#include "evo_ops.h"

using evo::CandidateSolution;
using evo::Matrix;

namespace {

int tournament_by_rank(const std::vector<int>& rank, const std::vector<double>& crowd, std::mt19937& rng) {
    std::uniform_int_distribution<int> dist(0, static_cast<int>(rank.size()) - 1);
    int a = dist(rng), b = dist(rng);
    if (rank[a] < rank[b]) return a;
    if (rank[b] < rank[a]) return b;
    if (crowd[a] > crowd[b]) return a;
    if (crowd[b] > crowd[a]) return b;
    return a;
}

}  // namespace

NSGAII::NSGAII(const Model& model, int pop_size)
    : model_(model), pop_size_(pop_size) {}

std::vector<CandidateSolution> NSGAII::run(const TerminationCriteria& term, ProgressLogger* logger, bool use_repair) {
    std::mt19937 rng(std::random_device{}());
    std::vector<Matrix> prod_pop = evo::random_population(model_, pop_size_, rng, 1.5, use_repair);
    std::vector<CandidateSolution> population = evo::evaluate_population(model_, prod_pop);

    int generation = 0;
    int evaluations = static_cast<int>(population.size());
    auto start = std::chrono::steady_clock::now();

    std::vector<int> rank(population.size(), 0);
    std::vector<double> crowd(population.size(), 0.0);
    {
        auto fronts = evo::fast_non_dominated_sort(population);
        for (size_t f = 0; f < fronts.size(); ++f) {
            auto cd = evo::crowding_distance(population, fronts[f]);
            for (size_t k = 0; k < fronts[f].size(); ++k) {
                rank[fronts[f][k]] = static_cast<int>(f);
                crowd[fronts[f][k]] = cd[k];
            }
        }
    }

    while (true) {
        if (term.max_generations > 0 && generation >= term.max_generations) break;
        if (term.max_evaluations > 0 && evaluations >= term.max_evaluations) break;
        if (term.max_seconds > 0.0) {
            double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            if (elapsed >= term.max_seconds) break;
        }

        std::vector<CandidateSolution> offspring;
        std::vector<Matrix> prod_child;
        offspring.reserve(pop_size_);
        prod_child.reserve(pop_size_);
        std::uniform_real_distribution<> ur(0.0, 1.0);
        while (static_cast<int>(offspring.size()) < pop_size_) {
            int a = std::uniform_int_distribution<int>(0, pop_size_ - 1)(rng);
            int b = std::uniform_int_distribution<int>(0, pop_size_ - 1)(rng);
            Matrix child_prod = evo::crossover_sbx(prod_pop[a], prod_pop[b], rng, 20.0, evo::ops_nsgaii.pc);
            if (ur(rng) < evo::ops_nsgaii.p_toggle) evo::mutate_toggle_on_off(model_, child_prod, rng);
            if (use_repair && ur(rng) < evo::ops_nsgaii.p_repair) evo::repair_demand(model_, child_prod);
            offspring.push_back(evo::evaluate(model_, child_prod));
            prod_child.push_back(std::move(child_prod));
            evaluations++;
        }

        std::vector<CandidateSolution> combined = population;
        combined.insert(combined.end(), offspring.begin(), offspring.end());
        std::vector<Matrix> prod_combined = prod_pop;
        prod_combined.insert(prod_combined.end(), prod_child.begin(), prod_child.end());

        auto fronts_all = evo::fast_non_dominated_sort(combined);
        std::vector<CandidateSolution> next;
        std::vector<Matrix> prod_next;
        std::vector<int> rank_next;
        std::vector<double> crowd_next;
        next.reserve(pop_size_);
        prod_next.reserve(pop_size_);
        rank_next.reserve(pop_size_);
        crowd_next.reserve(pop_size_);
        for (size_t f = 0; f < fronts_all.size() && static_cast<int>(next.size()) < pop_size_; ++f) {
            auto cd = evo::crowding_distance(combined, fronts_all[f]);
            std::vector<int> order(fronts_all[f].size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](int a, int b) { return cd[a] > cd[b]; });
            for (int idx_local : order) {
                int idx = fronts_all[f][idx_local];
                if (static_cast<int>(next.size()) >= pop_size_) break;
                next.push_back(combined[idx]);
                prod_next.push_back(prod_combined[idx]);
                rank_next.push_back(static_cast<int>(f));
                crowd_next.push_back(cd[idx_local]);
            }
        }
        population.swap(next);
        prod_pop.swap(prod_next);
        rank.swap(rank_next);
        crowd.swap(crowd_next);
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

    return population;
}
