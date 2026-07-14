#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "evo_ops.h"
#include "termination.h"

struct ExperimentConfig {
    std::vector<std::filesystem::path> instances;
    int repeats{1};
    int pop_size{60};
    TerminationCriteria term{};
    std::filesystem::path output_root{"out_results"};
    bool use_repair{true};
    evo::OperatorProbs ops_nsgaii{};
    evo::OperatorProbs ops_nsgaiii{};
    evo::OperatorProbs ops_moead{};
    evo::OperatorProbs ops_clsp{};
};

void run_experiments(const ExperimentConfig& cfg);
