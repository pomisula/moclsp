#include "experiment_runner.h"

#include <chrono>
#include <filesystem>
#include <iostream>
#include <string>

#include "algo_clsp.h"
#include "instance_loader.h"
#include "moead.h"
#include "nsgaii.h"
#include "nsgaiii.h"
#include "progress_logger.h"
#include "evo_ops.h"

namespace fs = std::filesystem;

void run_experiments(const ExperimentConfig& cfg) {
    for (const auto& inst_rel : cfg.instances) {
        try {
            Model model = load_instance_from_flat(fs::path(INSTANCE_ROOT) / inst_rel);
            evo::ops_nsgaii = cfg.ops_nsgaii;
            evo::ops_nsgaiii = cfg.ops_nsgaiii;
            evo::ops_moead = cfg.ops_moead;
            evo::ops_clsp = cfg.ops_clsp;
            TerminationCriteria term = cfg.term;
            if (term.max_evaluations <= 0) {
                term.max_evaluations = model.J * model.T * 5 * cfg.pop_size;
            }
            term.max_generations = 0;
            term.max_seconds = 0;
            fs::path run_root = cfg.output_root / inst_rel.string();

            for (int r = 0; r < cfg.repeats; ++r) {
                int run_id = r + 1;
                auto run_algo = [&](const std::string& name,
                                    auto runner) {
                    fs::path base = run_root / name / ("run" + std::to_string(run_id));
                    ProgressLogger logger("", "", "");
                    auto t0 = std::chrono::steady_clock::now();
                    auto pop = runner(logger);
                    double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

                    fs::path pop_file = base / (inst_rel.filename().string() + "_final_pop.csv");
                    logger.set_final_population(model, pop, pop_file.string());
                    (void)elapsed;
                };

                run_algo("NSGAII", [&](ProgressLogger& logger) {
                    NSGAII algo(model, cfg.pop_size);
                    return algo.run(term, &logger, cfg.use_repair);
                });
                run_algo("NSGAIII", [&](ProgressLogger& logger) {
                    NSGAIII algo(model, cfg.pop_size);
                    return algo.run(term, &logger, cfg.use_repair);
                });

                run_algo("MOEAD", [&](ProgressLogger& logger) {
                    MOEAD algo(model, cfg.pop_size);
                    return algo.run(term, &logger, cfg.use_repair);
                });
                run_algo("AlgoCLSP", [&](ProgressLogger& logger) {
                    AlgoCLSP algo(model, cfg.pop_size);
                    return algo.run(term, &logger);
                });
            }
        } catch (const std::exception& ex) {
            std::cerr << "Experiment skipped for " << inst_rel << ": " << ex.what() << "\n";
        }
    }
}
