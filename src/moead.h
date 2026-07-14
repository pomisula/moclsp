#pragma once

#include <vector>

#include "evo_ops.h"
#include "model.h"
#include "termination.h"
#include "progress_logger.h"

class MOEAD {
public:
    MOEAD(const Model& model, int pop_size = 60);
    std::vector<evo::CandidateSolution> run(const TerminationCriteria& term, ProgressLogger* logger = nullptr, bool use_repair = true);

private:
    const Model& model_;
    int pop_size_;
};
