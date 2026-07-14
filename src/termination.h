#pragma once

struct TerminationCriteria {
    int max_generations{200};
    int max_evaluations{0};
    double max_seconds{0.0};
};
