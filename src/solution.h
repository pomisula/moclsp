#pragma once

#include <vector>

#include "model.h"

class Solution {
public:
    explicit Solution(const Model& model)
        : production(model.J, std::vector<double>(model.T, 0.0)) {}

    void add_production(int j, int t, int /*m*/, double qty, double /*process_time*/, double /*setup_time*/) {
        if (qty <= 0) return;
        production[j][t] += qty;
    }

    const std::vector<std::vector<double>>& get_production() const { return production; }

private:
    std::vector<std::vector<double>> production;
};
