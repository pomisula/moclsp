# MOCLSP Test Instances

## About the Test Instances

This project contains 50 modified test instances based on classic CLSP benchmarks, designed for algorithm performance evaluation and comparative analysis.

### Instance Source

These instances are modified from the classic CLSP test set proposed by Trigeiro et al. (1989). The benchmark set is available at: https://suerie.de/testsets.html

**Reference:**
> Trigeiro, W. W., Thomas, L. J., & McClain, J. O. (1989). Capacitated Lot Sizing with Setup Times. *Management Science*, 35(3), 353-366. https://doi.org/10.1287/mnsc.35.3.353

### Instance Characteristics

The test set covers various problem scales with 5 instances for each configuration:

| Periods (T) | Items (J) | Number of Instances |
|------------|-----------|---------------------|
| 15         | 6         | 5 |
| 15         | 8         | 5 |
| 15         | 12        | 5 |
| 15         | 24        | 5 |
| 20         | 10        | 5 |
| 20         | 20        | 5 |
| 20         | 30        | 5 |
| 30         | 6         | 5 |
| 30         | 12        | 5 |
| 30         | 24        | 5 |

### Instance Data

Each instance contains the following key data:

- **INDEX**: Number of items (J), periods (T), and machines (M)
- **LAGKOST**: Inventory holding costs
- **L0/LT**: Initial/target inventory levels
- **KAPAZ**: Machine capacity constraints
- **P-BEDARF**: Product demands
- **PRODKOEF**: Production coefficients (processing time per unit)
- **RUESTK**: Setup costs
- **RUESTZ**: Setup times
- **UEBER-KS**: Overtime costs

### Naming Convention

Instance files follow the naming pattern:

```
inst_T{periods}_J{items}_{index}.txt
```

For example: `inst_T15_J12_13.txt` represents an instance with 15 periods and 12 items.

## Source Code

The source code and implementation details will be released publicly in the future.

---

*Last updated: January 2026*
