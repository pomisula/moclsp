#include <algorithm>
#include <filesystem>
#include <iostream>
#include <vector>

#include "experiment_runner.h"
#include "termination.h"

namespace fs = std::filesystem;

#ifndef INSTANCE_ROOT
#define INSTANCE_ROOT "."
#endif

namespace {

std::vector<fs::path> collect_instances(const fs::path& dir) {
    std::vector<fs::path> instances;
    if (!fs::exists(dir)) return instances;

    for (const auto& entry : fs::directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;
        if (entry.path().extension() != ".txt") continue;
        if (entry.path().filename().string().rfind("inst_", 0) != 0) continue;
        instances.push_back(fs::path("Instance") / entry.path().filename());
    }

    std::sort(instances.begin(), instances.end(),
              [](const fs::path& a, const fs::path& b) { return a.filename().string() < b.filename().string(); });
    return instances;
}

}  // namespace

int main() {
    try {
        const auto instances = collect_instances(fs::path(INSTANCE_ROOT) / "Instance");
        if (instances.empty()) {
            std::cerr << "No instance files found under: " << (fs::path(INSTANCE_ROOT) / "Instance") << "\n";
            return 1;
        }

        ExperimentConfig cfg;
        cfg.instances = instances;
        cfg.repeats = 20;
        cfg.pop_size = 100;
        cfg.output_root = fs::path(INSTANCE_ROOT) / "out_results";

        evo::OperatorProbs nsgaii_ops{};
        nsgaii_ops.pc = 0.5;
        nsgaii_ops.p_toggle = 0.7;
        nsgaii_ops.p_repair = 0.7;

        evo::OperatorProbs nsgaiii_ops{};
        nsgaiii_ops.pc = 0.5;
        nsgaiii_ops.p_toggle = 0.3;
        nsgaiii_ops.p_repair = 1.0;

        evo::OperatorProbs moead_ops{};
        moead_ops.pc = 0.3;
        moead_ops.p_toggle = 0.3;
        moead_ops.p_repair = 0.7;

        evo::OperatorProbs clsp_ops{};
        clsp_ops.pc = 1.0;
        clsp_ops.p_shift = 1.0;
        clsp_ops.p_toggle = 0.7;
        clsp_ops.p_repair = 0.7;
        clsp_ops.lns_iters = 1;
        clsp_ops.lns_relax_frac = 0.1;
        clsp_ops.lns_clusters = 4;

        cfg.ops_nsgaii = nsgaii_ops;
        cfg.ops_nsgaiii = nsgaiii_ops;
        cfg.ops_moead = moead_ops;
        cfg.ops_clsp = clsp_ops;

        TerminationCriteria term;
        term.max_generations = 0;
        term.max_evaluations = 0;
        term.max_seconds = 0.0;
        cfg.term = term;

        run_experiments(cfg);
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
