#include "evo_ops.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace evo {

OperatorProbs ops_nsgaii{};
OperatorProbs ops_nsgaiii{};
OperatorProbs ops_moead{};
OperatorProbs ops_clsp{};

void repair_demand(const Model& model, Matrix& prod) {
    std::vector<double> inv(model.J, 0.0);
    if (static_cast<int>(model.initial_inventory.size()) == model.J) {
        inv = model.initial_inventory;
    }
    for (int j = 0; j < model.J; ++j) {
        for (int t = 0; t < model.T; ++t) {
            double dem = model.demand[j][t];
            double avail = inv[j] + prod[j][t];
            if (avail >= dem) {
                inv[j] = avail - dem;
            } else {
                double need = dem - avail;
                prod[j][t] += need;
                inv[j] = 0.0;
            }
        }
    }
}

std::vector<Matrix> random_population(const Model& model, int pop_size, std::mt19937& rng, double scale,
                                      bool use_repair) {
    std::vector<Matrix> pop;
    pop.reserve(pop_size);
    for (int i = 0; i < pop_size; ++i) {
        Matrix p(model.J, std::vector<double>(model.T, 0.0));
        for (int j = 0; j < model.J; ++j) {
            for (int t = 0; t < model.T; ++t) {
                double cap = std::max(1.0, model.demand[j][t] * scale);
                p[j][t] = std::uniform_real_distribution<double>(0.0, cap)(rng);
            }
        }
        if (use_repair) repair_demand(model, p);
        pop.push_back(std::move(p));
    }
    return pop;
}

Solution build_solution(const Model& model, const Matrix& prod) {
    Solution sol(model);
    for (int j = 0; j < model.J; ++j) {
        for (int t = 0; t < model.T; ++t) {
            double qty = prod[j][t];
            if (qty > 0.0) {
                sol.add_production(j, t, 0, qty, model.prod_coeff[0][j], model.setup_time[0][j]);
            }
        }
    }
    return sol;
}

CandidateSolution evaluate(const Model& model, const Matrix& prod) {
    CandidateSolution res(model);
    res.plan = build_solution(model, prod);
    res.metrics = model.evaluate_solution(res.plan);
    return res;
}

std::vector<CandidateSolution> evaluate_population(const Model& model, const std::vector<Matrix>& prod_vec) {
    std::vector<CandidateSolution> res;
    res.reserve(prod_vec.size());
    for (const auto& p : prod_vec) {
        res.push_back(evaluate(model, p));
    }
    return res;
}

Matrix crossover_sbx(const Matrix& a, const Matrix& b, std::mt19937& rng, double eta_c, double pc) {
    Matrix c = a;
    std::uniform_real_distribution<> ur(0.0, 1.0);
    if (ur(rng) > pc) return c;
    for (size_t j = 0; j < a.size(); ++j) {
        for (size_t t = 0; t < a[j].size(); ++t) {
            double u = ur(rng);
            double beta;
            if (u <= 0.5) {
                beta = std::pow(2 * u, 1.0 / (eta_c + 1));
            } else {
                beta = std::pow(1.0 / (2.0 * (1.0 - u)), 1.0 / (eta_c + 1));
            }
            c[j][t] = 0.5 * ((1 + beta) * a[j][t] + (1 - beta) * b[j][t]);
        }
    }
    return c;
}

void mutate_shift(const Model& model, Matrix& prod, std::mt19937& rng) {
    if (model.T <= 1 || model.J <= 0) return;
    std::uniform_int_distribution<int> item_dist(0, model.J - 1);
    int j = item_dist(rng);
    std::uniform_int_distribution<int> t_dist(1, model.T - 1);
    int t = t_dist(rng);
    double qty = prod[j][t];
    if (qty <= 0.0) return;
    double move = std::max(1.0, qty / 2.0);
    std::uniform_int_distribution<int> s_dist(0, t - 1);
    int s = s_dist(rng);
    prod[j][t] -= move;
    prod[j][s] += move;
}

void mutate_poly(Matrix& prod, std::mt19937& rng, double eta_m, double pm) {
    int J = static_cast<int>(prod.size());
    if (J == 0) return;
    int T = static_cast<int>(prod[0].size());
    if (pm < 0) pm = 1.0 / static_cast<double>(J * T);
    std::uniform_real_distribution<> ur(0.0, 1.0);
    for (int j = 0; j < J; ++j) {
        for (int t = 0; t < T; ++t) {
            if (ur(rng) > pm) continue;
            double y = prod[j][t];
            double yl = 0.0;
            double yu = std::max(1.0, y + 10.0);
            double delta1 = (y - yl) / (yu - yl);
            double delta2 = (yu - y) / (yu - yl);
            double rnd = ur(rng);
            double mut_pow = 1.0 / (eta_m + 1.0);
            double deltaq;
            if (rnd <= 0.5) {
                double xy = 1.0 - delta1;
                double val = 2.0 * rnd + (1.0 - 2.0 * rnd) * std::pow(xy, eta_m + 1.0);
                deltaq = std::pow(val, mut_pow) - 1.0;
            } else {
                double xy = 1.0 - delta2;
                double val = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * std::pow(xy, eta_m + 1.0);
                deltaq = 1.0 - std::pow(val, mut_pow);
            }
            y = y + deltaq * (yu - yl);
            y = std::min(std::max(y, yl), yu);
            prod[j][t] = y;
        }
    }
}

void mutate_toggle_on_off(const Model& model, Matrix& prod, std::mt19937& rng,
                          double p_on, double p_off, double scale) {
    int dim = model.J * model.T;
    double p_on_use = (p_on < 0) ? (1.0 / std::max(dim, 1)) : p_on;
    double p_off_use = (p_off < 0) ? (1.0 / std::max(dim, 1)) : p_off;
    std::uniform_real_distribution<> ur(0.0, 1.0);
    for (int j = 0; j < model.J; ++j) {
        for (int t = 0; t < model.T; ++t) {
            double u = ur(rng);
            if (prod[j][t] <= 0.0 && u < p_on_use) {
                double cap = std::max(1.0, model.demand[j][t] * scale);
                prod[j][t] = std::uniform_real_distribution<double>(0.0, cap)(rng);
            } else if (prod[j][t] > 0.0 && u < p_off_use) {
                prod[j][t] = 0.0;
            }
        }
    }
}

bool dominates(const SolutionMetrics& a, const SolutionMetrics& b, double tol) {
    bool feas_a = a.unmet_demand <= tol;
    bool feas_b = b.unmet_demand <= tol;
    if (feas_a && !feas_b) return true;
    if (!feas_a && feas_b) return false;
    if (!feas_a && !feas_b) return a.unmet_demand < b.unmet_demand - tol;
    bool better_or_eq = (a.holding_cost <= b.holding_cost + tol) &&
                        (a.setup_cost <= b.setup_cost + tol) &&
                        (a.overtime_cost <= b.overtime_cost + tol) &&
                        (a.unmet_demand <= b.unmet_demand + tol);
    bool strict = (a.holding_cost < b.holding_cost - tol) ||
                  (a.setup_cost < b.setup_cost - tol) ||
                  (a.overtime_cost < b.overtime_cost - tol) ||
                  (a.unmet_demand < b.unmet_demand - tol);
    return better_or_eq && strict;
}

std::vector<int> non_dominated_indices(const std::vector<SolutionMetrics>& ms, double tol) {
    std::vector<int> nd;
    for (int i = 0; i < static_cast<int>(ms.size()); ++i) {
        bool dom = false;
        for (int j = 0; j < static_cast<int>(ms.size()); ++j) {
            if (i == j) continue;
            if (dominates(ms[j], ms[i], tol)) {
                dom = true;
                break;
            }
        }
        if (!dom) nd.push_back(i);
    }
    return nd;
}

}  // namespace evo
