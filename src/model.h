#pragma once

#include <string>
#include <vector>

class Solution;

struct SolutionMetrics {
    double total_cost{0.0};
    double holding_cost{0.0};
    double setup_cost{0.0};
    double overtime_cost{0.0};
    double max_capacity_overuse{0.0};  // after applying overtime this should stay ~0
    double unmet_demand{0.0};
};

struct Model {
    int J{};  // items
    int T{};  // periods
    int M{};  // machines

    std::vector<std::vector<double>> demand;       // J x T
    std::vector<std::vector<double>> capacity;     // M x T
    std::vector<double> holding_cost;              // J
    std::vector<double> setup_cost;                // J
    std::vector<double> initial_inventory;         // J
    std::vector<std::vector<double>> prod_coeff;   // M x J (time per unit)
    std::vector<std::vector<double>> setup_time;   // M x J
    std::vector<double> overtime_cost;             // M (optional, may be empty)

    const SolutionMetrics& evaluate_solution(const Solution& sol) const;

private:
    mutable SolutionMetrics last_metrics_{};
};
