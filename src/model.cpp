#include "model.h"

#include <algorithm>

#include "solution.h"

const SolutionMetrics& Model::evaluate_solution(const Solution& sol) const {
    SolutionMetrics metrics;

    const auto& prod = sol.get_production();

    std::vector<double> used_time(T, 0.0);
    for (int t = 0; t < T; ++t) {
        for (int j = 0; j < J; ++j) {
            double qty = prod[j][t];
            if (qty <= 0.0) continue;
            used_time[t] += qty * prod_coeff[0][j];
            used_time[t] += setup_time[0][j];
            metrics.setup_cost += setup_cost[j];
        }
    }

    for (int t = 0; t < T; ++t) {
        double over = used_time[t] - capacity[0][t];
        if (over > 0.0) {
            metrics.max_capacity_overuse = std::max(metrics.max_capacity_overuse, over);
            if (!overtime_cost.empty()) metrics.overtime_cost += over * overtime_cost[0];
        }
    }

    std::vector<double> inv = initial_inventory;
    if (static_cast<int>(inv.size()) != J) inv.assign(J, 0.0);
    for (int t = 0; t < T; ++t) {
        for (int j = 0; j < J; ++j) {
            double avail = inv[j] + prod[j][t];
            double dem = demand[j][t];
            if (avail >= dem) {
                inv[j] = avail - dem;
            } else {
                double unmet = dem - avail;
                inv[j] = 0.0;
                metrics.unmet_demand += unmet;
            }
            metrics.holding_cost += inv[j] * holding_cost[j];
        }
    }

    metrics.total_cost = metrics.holding_cost + metrics.setup_cost + metrics.overtime_cost;
    last_metrics_ = metrics;
    return last_metrics_;
}
