# MOCLSP

## About the Test Instances

MOCLSP (**Multi-Objective Capacitated Lot Sizing Problem**).

This project contains 50 modified test instances based on classic CLSP benchmarks, designed for algorithm performance evaluation and comparative analysis.

The MOCLSP instances used in this project are provided in `Instance.7z`.

### Instance Source

These instances are derived from the classic CLSP test set proposed by Trigeiro et al. (1989). The original benchmark set is available at: https://suerie.de/testsets.html

The instance set in this repository follows the original benchmark structure and uses the same problem parameters and definitions as the Trigeiro CLSP instances, including demand, capacity, holding cost, production coefficients, setup cost, setup time, and initial inventory information. The selected instances cover the period-item scale combinations used in the original benchmark data that are relevant to this study.

**Reference:**
> Trigeiro, W. W., Thomas, L. J., & McClain, J. O. (1989). Capacitated Lot Sizing with Setup Times. *Management Science*, 35(3), 353-366. https://doi.org/10.1287/mnsc.35.3.353

### Instance Characteristics

The original benchmark data contains 10 relevant period-item (`T x J`) configurations. This test set covers all 10 configurations. For each configuration, 5 corresponding instances are included because 5 is the smallest number of available instances among these configurations in the original benchmark data. This keeps the instance set balanced across problem scales:

| Periods (T) | Items (J) | Number of Instances |
| ----------- | --------- | ------------------- |
| 15          | 6         | 5                   |
| 15          | 8         | 5                   |
| 15          | 12        | 5                   |
| 15          | 24        | 5                   |
| 20          | 10        | 5                   |
| 20          | 20        | 5                   |
| 20          | 30        | 5                   |
| 30          | 6         | 5                   |
| 30          | 12        | 5                   |
| 30          | 24        | 5                   |

### Instance Data Sections

Each instance file is organized into labeled data sections. The labels follow the naming convention of the original benchmark data:

| Section | Description | Structure |
| ------- | ----------- | --------- |
| `INDEX` | Instance dimensions, number of items (`J`), periods (`T`), and machines/resources (`M`) | 1 row with 3 values, `J, T, M` |
| `LAGKOST` | Unit holding cost for each item | 1 row with `J` values |
| `L0` | Beginning inventory for each item | 1 row with `J` values |
| `LT` | Ending inventory for each item | 1 row with `J` values |
| `KAPAZ` | Available capacity for each machine/resource and period | `M` rows and `T` columns |
| `P-BEDARF` | Demand for each item and period | `J` rows and `T` columns |
| `PRODKOEF` | Production coefficient, capacity needed on a machine/resource for one unit of an item | One row per machine/resource-item pair, `m, j, value` |
| `RUESTK` | Setup cost for each item | `J` rows and 1 column |
| `RUESTZ` | Setup time for each machine/resource and item | One row per machine/resource-item pair, `m, j, value` |
| `UEBER-KS` | Overtime cost per unit of capacity overuse | `M` rows and 1 column |
| `DIREKT-B` | Demand coefficients preserved from the original instance format | Empty for the included single-level instances |

### Naming Convention

Instance files follow the naming pattern:

```
inst_T{periods}_J{items}_{index}.txt
```

For example: `inst_T15_J12_13.txt` represents an instance with 15 periods and 12 items.

## Code and Run

Implemented algorithms:

- `AlgoCLSP`
- `NSGAII`
- `NSGAIII`
- `MOEAD`

Environment:

- CMake >= 3.15
- C++17 compiler

Main experiment settings are mainly in `src/main.cpp` and `src/experiment_runner.cpp`.
