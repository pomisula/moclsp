#pragma once

#include <filesystem>
#include <fstream>
#include <string>

#include "model.h"

class ProgressLogger {
public:
    ProgressLogger(const std::string& metrics_path,
                   const std::string& pop_dir = "",
                   const std::string& pop_metrics_path_or_dir = "") {
        (void)metrics_path;
        (void)pop_dir;
        (void)pop_metrics_path_or_dir;
    }

    ~ProgressLogger() = default;

    void log(int generation, int evals, double time_sec, const SolutionMetrics& m) {
        last_generation_ = generation;
        last_evaluations_ = evals;
        last_time_ = time_sec;
        last_metrics_ = m;
    }

    template <typename Vec>
    void log_population_metrics(int generation, const Vec& pop) {
        (void)generation;
        (void)pop;
    }

    template <typename Vec>
    void log_population(int generation, const Model& model, const Vec& pop) {
        (void)generation;
        (void)model;
        (void)pop;
    }

    template <typename Vec>
    void set_final_population(const Model& model, const Vec& pop, const std::string& path) {
        std::filesystem::path out(path);
        ensure_parent(out);
        std::ofstream ofs(out, std::ios::trunc);
        if (!ofs.is_open()) return;

        ofs << "evaluations,time_sec,idx,total,hold,setup,overtime,unmet";
        for (int j = 0; j < model.J; ++j) {
            for (int t = 0; t < model.T; ++t) ofs << ",x_" << j << "_" << t;
        }
        ofs << "\n";

        for (size_t i = 0; i < pop.size(); ++i) {
            const auto& m = pop[i].metrics;
            const auto& prod = pop[i].plan.get_production();
            ofs << last_evaluations_ << "," << last_time_ << "," << i << "," << m.total_cost << ","
                << m.holding_cost << "," << m.setup_cost << "," << m.overtime_cost << "," << m.unmet_demand;
            for (int j = 0; j < model.J; ++j) {
                for (int t = 0; t < model.T; ++t) ofs << "," << prod[j][t];
            }
            ofs << "\n";
        }
    }

    void flush() {}

    int last_generation() const { return last_generation_; }
    int last_evaluations() const { return last_evaluations_; }
    double last_time() const { return last_time_; }
    const SolutionMetrics& last_metrics() const { return last_metrics_; }

private:
    static void ensure_parent(const std::filesystem::path& p) {
        auto parent = p.parent_path();
        if (!parent.empty()) {
            std::filesystem::create_directories(parent);
        }
    }

    int last_generation_{0};
    int last_evaluations_{0};
    double last_time_{0.0};
    SolutionMetrics last_metrics_{};
};
