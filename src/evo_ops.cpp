#include "evo_ops.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace evo {

OperatorProbs ops_nsgaii{1.0, 0.6, 0.3, 1.0, 0.0, 3, 0.25, 0};
OperatorProbs ops_nsgaiii{0.5, 0.0, 0.3, 1.0, 0.0, 3, 0.25, 0};
OperatorProbs ops_moead{1.0, 0.0, 0.3, 1.0, 0.0, 3, 0.25, 0};
OperatorProbs ops_clsp{1.0, 0.6, 0.3, 1.0, 0.5, 3, 0.25, 0};

namespace {

int comb2(int n) {
    return n * (n - 1) / 2;
}

std::vector<std::array<double, 3>> simplex_lattice_3d(int h) {
    std::vector<std::array<double, 3>> points;
    if (h <= 0) return points;
    points.reserve(comb2(h + 2));
    for (int a = 0; a <= h; ++a) {
        for (int b = 0; b <= h - a; ++b) {
            int c = h - a - b;
            double s = static_cast<double>(h);
            points.push_back({a / s, b / s, c / s});
        }
    }
    return points;
}

}  // namespace

std::vector<std::array<double, 3>> uniform_points_nbi_3d(int n) {
    constexpr int m = 3;
    int h1 = 1;
    while (comb2(h1 + m) <= n) ++h1;

    std::vector<std::array<double, 3>> points = simplex_lattice_3d(h1);
    if (h1 < m) {
        int h2 = 0;
        while (static_cast<int>(points.size()) + comb2(h2 + m) <= n) ++h2;
        if (h2 > 0) {
            auto second = simplex_lattice_3d(h2);
            points.reserve(points.size() + second.size());
            for (auto p : second) {
                for (double& v : p) v = v / 2.0 + 1.0 / (2.0 * m);
                points.push_back(p);
            }
        }
    }

    for (auto& p : points) {
        for (double& v : p) v = std::max(v, 1e-6);
    }
    return points;
}

Matrix demand_matrix(const Model& model) {
    Matrix dm(model.J, std::vector<double>(model.T, 0.0));
    for (int j = 0; j < model.J; ++j) {
        for (int t = 0; t < model.T; ++t) {
            dm[j][t] = model.demand[j][t];
        }
    }
    return dm;
}

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

Matrix crossover_mix(const Matrix& a, const Matrix& b, std::mt19937& rng) {
    Matrix c = a;
    for (size_t j = 0; j < a.size(); ++j) {
        if (std::uniform_real_distribution<>(0.0, 1.0)(rng) < 0.5) {
            c[j] = b[j];
        }
    }
    return c;
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

Matrix crossover_uniform_cells(const Matrix& a, const Matrix& b, std::mt19937& rng, double p_b) {
    Matrix c = a;
    std::uniform_real_distribution<> ur(0.0, 1.0);
    for (size_t j = 0; j < a.size(); ++j) {
        for (size_t t = 0; t < a[j].size(); ++t) {
            if (ur(rng) < p_b) c[j][t] = b[j][t];
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

void mutate_zero_sparse(const Model& model, Matrix& prod, std::mt19937& rng, double drop_prob) {
    int dim = model.J * model.T;
    double p = (drop_prob < 0) ? (1.0 / std::max(dim, 1)) : drop_prob;
    std::uniform_real_distribution<> ur(0.0, 1.0);
    for (int j = 0; j < model.J; ++j) {
        for (int t = 0; t < model.T; ++t) {
            if (ur(rng) < p) prod[j][t] = 0.0;
        }
    }
}

void mutate_resample(const Model& model, Matrix& prod, std::mt19937& rng, double resample_prob, double scale) {
    int dim = model.J * model.T;
    double p = (resample_prob < 0) ? (1.0 / std::max(dim, 1)) : resample_prob;
    std::uniform_real_distribution<> ur(0.0, 1.0);
    for (int j = 0; j < model.J; ++j) {
        for (int t = 0; t < model.T; ++t) {
            if (ur(rng) < p) {
                double cap = std::max(1.0, model.demand[j][t] * scale);
                prod[j][t] = std::uniform_real_distribution<double>(0.0, cap)(rng);
            }
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

void lns_relax(const Model& model, Matrix& prod, std::mt19937& rng, int relax_periods) {
    if (relax_periods < 0) relax_periods = std::max(1, model.T / 4);
    std::uniform_int_distribution<int> t_dist(0, model.T - 1);
    std::vector<int> chosen;
    chosen.reserve(relax_periods);
    while (static_cast<int>(chosen.size()) < relax_periods) {
        int t = t_dist(rng);
        if (std::find(chosen.begin(), chosen.end(), t) == chosen.end()) chosen.push_back(t);
    }
    for (int t : chosen) {
        for (int j = 0; j < model.J; ++j) prod[j][t] = 0;
    }
    repair_demand(model, prod);
}

int tournament_select(const std::vector<double>& score, int k, std::mt19937& rng) {
    std::uniform_int_distribution<int> dist(0, static_cast<int>(score.size()) - 1);
    int best = dist(rng);
    for (int i = 1; i < k; ++i) {
        int cand = dist(rng);
        if (score[cand] < score[best]) best = cand;
    }
    return best;
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
