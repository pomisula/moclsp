#include "nsgaiii.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

using evo::CandidateSolution;
using evo::Matrix;

namespace {

constexpr double kTol = 1e-9;
constexpr double kEps = 1e-12;

std::array<double, 3> obj3(const CandidateSolution& s) {
    return {s.metrics.holding_cost, s.metrics.setup_cost, s.metrics.overtime_cost};
}

double l2norm3(const std::array<double, 3>& v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

double dot3(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

std::vector<std::array<double, 3>> uniform_points_3d(int n) {
    int h = 1;
    while ((h + 1) * (h + 2) / 2 < n) ++h;
    std::vector<std::array<double, 3>> z;
    for (int a = 0; a <= h; ++a) {
        for (int b = 0; b <= h - a; ++b) {
            int c = h - a - b;
            double s = static_cast<double>(h);
            z.push_back({a / s, b / s, c / s});
        }
    }
    if (static_cast<int>(z.size()) > n) z.resize(n);
    while (static_cast<int>(z.size()) < n) {
        double r1 = static_cast<double>(rand()) / RAND_MAX;
        double r2 = static_cast<double>(rand()) / RAND_MAX;
        double r3 = static_cast<double>(rand()) / RAND_MAX;
        double sum = r1 + r2 + r3 + kEps;
        z.push_back({r1 / sum, r2 / sum, r3 / sum});
    }
    for (auto& v : z) {
        for (double& x : v) x = std::max(x, 1e-6);
    }
    return z;
}

bool solve3x3(const std::array<std::array<double, 3>, 3>& a, const std::array<double, 3>& b,
              std::array<double, 3>& x) {
    double m[3][4] = {
        {a[0][0], a[0][1], a[0][2], b[0]},
        {a[1][0], a[1][1], a[1][2], b[1]},
        {a[2][0], a[2][1], a[2][2], b[2]},
    };
    for (int col = 0; col < 3; ++col) {
        int piv = col;
        for (int r = col + 1; r < 3; ++r) {
            if (std::abs(m[r][col]) > std::abs(m[piv][col])) piv = r;
        }
        if (std::abs(m[piv][col]) < kEps) return false;
        if (piv != col) {
            for (int c = col; c < 4; ++c) std::swap(m[col][c], m[piv][c]);
        }
        double d = m[col][col];
        for (int c = col; c < 4; ++c) m[col][c] /= d;
        for (int r = 0; r < 3; ++r) {
            if (r == col) continue;
            double f = m[r][col];
            for (int c = col; c < 4; ++c) m[r][c] -= f * m[col][c];
        }
    }
    x = {m[0][3], m[1][3], m[2][3]};
    return true;
}

std::vector<int> last_selection(const std::vector<std::array<double, 3>>& pop_obj1,
                                const std::vector<std::array<double, 3>>& pop_obj2,
                                int k,
                                const std::vector<std::array<double, 3>>& z,
                                const std::array<double, 3>& zmin,
                                std::mt19937& rng) {
    std::vector<int> chosen(pop_obj2.size(), 0);
    if (k <= 0 || pop_obj2.empty()) return chosen;

    std::vector<std::array<double, 3>> pop_obj;
    pop_obj.reserve(pop_obj1.size() + pop_obj2.size());
    for (auto p : pop_obj1) pop_obj.push_back({p[0] - zmin[0], p[1] - zmin[1], p[2] - zmin[2]});
    for (auto p : pop_obj2) pop_obj.push_back({p[0] - zmin[0], p[1] - zmin[1], p[2] - zmin[2]});

    int n = static_cast<int>(pop_obj.size());
    int n1 = static_cast<int>(pop_obj1.size());
    int nz = static_cast<int>(z.size());

    std::array<int, 3> extreme = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {
        std::array<double, 3> w = {1e-6, 1e-6, 1e-6};
        w[i] += 1.0;
        double best = std::numeric_limits<double>::infinity();
        int best_idx = 0;
        for (int p = 0; p < n; ++p) {
            double asf = std::max({pop_obj[p][0] / w[0], pop_obj[p][1] / w[1], pop_obj[p][2] / w[2]});
            if (asf < best) {
                best = asf;
                best_idx = p;
            }
        }
        extreme[i] = best_idx;
    }

    std::array<double, 3> a = {0.0, 0.0, 0.0};
    {
        std::array<std::array<double, 3>, 3> h = {
            pop_obj[extreme[0]], pop_obj[extreme[1]], pop_obj[extreme[2]]
        };
        std::array<double, 3> one = {1.0, 1.0, 1.0};
        std::array<double, 3> hp{};
        bool ok = solve3x3(h, one, hp);
        if (ok) {
            for (int d = 0; d < 3; ++d) {
                if (std::abs(hp[d]) < kEps) ok = false;
            }
            if (ok) {
                for (int d = 0; d < 3; ++d) a[d] = 1.0 / hp[d];
            }
        }
        if (!ok) {
            a = {0.0, 0.0, 0.0};
            for (const auto& p : pop_obj) {
                for (int d = 0; d < 3; ++d) a[d] = std::max(a[d], p[d]);
            }
        }
        for (int d = 0; d < 3; ++d) {
            if (a[d] <= kEps || std::isnan(a[d]) || std::isinf(a[d])) a[d] = 1.0;
        }
    }

    for (auto& p : pop_obj) {
        for (int d = 0; d < 3; ++d) p[d] /= a[d];
    }

    std::vector<int> pi(n, 0);
    std::vector<double> dist(n, std::numeric_limits<double>::infinity());
    for (int i = 0; i < n; ++i) {
        double np = l2norm3(pop_obj[i]);
        if (np < kEps) np = kEps;
        for (int j = 0; j < nz; ++j) {
            double nzj = l2norm3(z[j]);
            if (nzj < kEps) nzj = kEps;
            double cosine = dot3(pop_obj[i], z[j]) / (np * nzj);
            cosine = std::max(-1.0, std::min(1.0, cosine));
            double d = np * std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
            if (d < dist[i]) {
                dist[i] = d;
                pi[i] = j;
            }
        }
    }

    std::vector<int> rho(nz, 0);
    for (int i = 0; i < n1; ++i) rho[pi[i]]++;

    std::vector<int> z_choose(nz, 1);
    int picked = 0;
    while (picked < k) {
        std::vector<int> active;
        active.reserve(nz);
        for (int j = 0; j < nz; ++j) {
            if (z_choose[j]) active.push_back(j);
        }
        if (active.empty()) break;

        int min_rho = std::numeric_limits<int>::max();
        for (int j : active) min_rho = std::min(min_rho, rho[j]);
        std::vector<int> jmin;
        for (int j : active) {
            if (rho[j] == min_rho) jmin.push_back(j);
        }
        int j = jmin[std::uniform_int_distribution<int>(0, static_cast<int>(jmin.size()) - 1)(rng)];

        std::vector<int> i_set;
        for (int i = 0; i < static_cast<int>(pop_obj2.size()); ++i) {
            if (!chosen[i] && pi[n1 + i] == j) i_set.push_back(i);
        }
        if (!i_set.empty()) {
            int s = 0;
            if (rho[j] == 0) {
                double best_d = std::numeric_limits<double>::infinity();
                for (int r = 0; r < static_cast<int>(i_set.size()); ++r) {
                    double d = dist[n1 + i_set[r]];
                    if (d < best_d) {
                        best_d = d;
                        s = r;
                    }
                }
            } else {
                s = std::uniform_int_distribution<int>(0, static_cast<int>(i_set.size()) - 1)(rng);
            }
            chosen[i_set[s]] = 1;
            rho[j]++;
            picked++;
        } else {
            z_choose[j] = 0;
        }
    }
    return chosen;
}

std::vector<CandidateSolution> environmental_selection(const std::vector<CandidateSolution>& combined, int n_select,
                                                       const std::vector<std::array<double, 3>>& z,
                                                       const std::array<double, 3>& zmin,
                                                       std::mt19937& rng) {
    auto fronts = evo::fast_non_dominated_sort(combined);
    std::vector<CandidateSolution> next;
    next.reserve(n_select);
    if (fronts.empty()) return next;

    int max_f = -1;
    int count = 0;
    for (int f = 0; f < static_cast<int>(fronts.size()); ++f) {
        if (count + static_cast<int>(fronts[f].size()) <= n_select) {
            count += static_cast<int>(fronts[f].size());
            max_f = f;
        } else {
            max_f = f;
            break;
        }
    }

    std::vector<int> fixed_idx;
    std::vector<int> last_idx;
    if (count == n_select) {
        for (int f = 0; f <= max_f; ++f) {
            fixed_idx.insert(fixed_idx.end(), fronts[f].begin(), fronts[f].end());
        }
        for (int idx : fixed_idx) next.push_back(combined[idx]);
        return next;
    }

    for (int f = 0; f < max_f; ++f) {
        fixed_idx.insert(fixed_idx.end(), fronts[f].begin(), fronts[f].end());
    }
    if (max_f >= 0 && max_f < static_cast<int>(fronts.size())) {
        last_idx = fronts[max_f];
    }

    std::vector<std::array<double, 3>> pop1;
    std::vector<std::array<double, 3>> pop2;
    pop1.reserve(fixed_idx.size());
    pop2.reserve(last_idx.size());
    for (int idx : fixed_idx) pop1.push_back(obj3(combined[idx]));
    for (int idx : last_idx) pop2.push_back(obj3(combined[idx]));

    int k = n_select - static_cast<int>(fixed_idx.size());
    auto choose = last_selection(pop1, pop2, k, z, zmin, rng);

    for (int idx : fixed_idx) next.push_back(combined[idx]);
    for (int i = 0; i < static_cast<int>(last_idx.size()); ++i) {
        if (choose[i]) next.push_back(combined[last_idx[i]]);
    }
    return next;
}

std::array<double, 3> compute_zmin_feasible(const std::vector<CandidateSolution>& pop) {
    std::array<double, 3> zmin = {std::numeric_limits<double>::infinity(),
                                  std::numeric_limits<double>::infinity(),
                                  std::numeric_limits<double>::infinity()};
    bool any_feasible = false;
    for (const auto& s : pop) {
        if (s.metrics.unmet_demand <= kTol) {
            auto o = obj3(s);
            for (int d = 0; d < 3; ++d) zmin[d] = std::min(zmin[d], o[d]);
            any_feasible = true;
        }
    }
    if (!any_feasible) zmin = {1.0, 1.0, 1.0};
    return zmin;
}

}  // namespace

NSGAIII::NSGAIII(const Model& model, int pop_size) : model_(model), pop_size_(pop_size) {}

std::vector<CandidateSolution> NSGAIII::run(const TerminationCriteria& term, ProgressLogger* logger, bool use_repair) {
    std::mt19937 rng(std::random_device{}());

    auto z = uniform_points_3d(pop_size_);
    std::vector<Matrix> prod_pop = evo::random_population(model_, pop_size_, rng, 1.5, use_repair);
    std::vector<CandidateSolution> population = evo::evaluate_population(model_, prod_pop);
    auto zmin = compute_zmin_feasible(population);

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

        std::vector<CandidateSolution> offspring;
        std::vector<Matrix> prod_child;
        offspring.reserve(pop_size_);
        prod_child.reserve(pop_size_);
        std::uniform_real_distribution<> ur(0.0, 1.0);
        std::uniform_int_distribution<int> psel(0, pop_size_ - 1);
        while (static_cast<int>(offspring.size()) < pop_size_) {
            int a = psel(rng);
            int b = psel(rng);
            Matrix child_prod = evo::crossover_sbx(prod_pop[a], prod_pop[b], rng, 20.0, evo::ops_nsgaiii.pc);
            if (ur(rng) < evo::ops_nsgaiii.p_toggle) evo::mutate_toggle_on_off(model_, child_prod, rng);
            if (use_repair && ur(rng) < evo::ops_nsgaiii.p_repair) evo::repair_demand(model_, child_prod);
            offspring.push_back(evo::evaluate(model_, child_prod));
            prod_child.push_back(std::move(child_prod));
            evaluations++;
        }

        std::vector<CandidateSolution> combined = population;
        combined.insert(combined.end(), offspring.begin(), offspring.end());
        std::vector<Matrix> prod_combined = prod_pop;
        prod_combined.insert(prod_combined.end(), prod_child.begin(), prod_child.end());

        zmin = compute_zmin_feasible(combined);
        auto next = environmental_selection(combined, pop_size_, z, zmin, rng);

        std::vector<Matrix> next_prod;
        next_prod.reserve(next.size());
        for (const auto& cs : next) {
            next_prod.push_back(cs.plan.get_production());
        }
        population.swap(next);
        prod_pop.swap(next_prod);
        generation++;

        if (logger && !population.empty()) {
            const auto* best = &population.front().metrics;
            for (const auto& v : population) {
                if (v.metrics.unmet_demand < best->unmet_demand - kTol ||
                    (std::abs(v.metrics.unmet_demand - best->unmet_demand) < kTol &&
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
