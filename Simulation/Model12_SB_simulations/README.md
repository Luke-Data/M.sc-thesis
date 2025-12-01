# Model 1.2 Smoothed Backfitting Simulations

## Simulation Design (2×2 Factorial)

Monte Carlo study evaluating wsbackfit performance under different conditions.

### Fixed Elements

- **Model**: Y = X₁f₁(Z₁) + X₂f₂(Z₁) + X₃[f₃(Z₂) + f₄(Z₁)] + ε, where ε ~ N(0,1)
- **Predictors**:
  - X₁ ~ N(2, 3)
  - X₂ ~ Exp(3)
  - X₃ ~ Gamma(1, 1)
- **True coefficient functions**:
  - f₁(z) = 2.5 sin(1.5z) exp(-0.2z²)
  - f₂(z) = z cos(2z)
  - f₃(z) = sin(1.5z)
  - f₄(z) = 3/(1+z²) cos(2z)
- **Monte Carlo replications**: nsim = 100
- **Estimation method**: wsbackfit with automatic bandwidth selection (h = -1)
- **Evaluation grid**: x_grid = seq(-3, 3, length = 500)

### Experimental Factors

| Factor | Levels | Description |
|--------|--------|-------------|
| **Z Distribution** | Uniform, Normal | Smoothing variables Z₁, Z₂ |
| **Sample Size** | 500, 1000 | Number of observations |

### Four Scenarios

| Scenario | Z₁, Z₂ Distribution | Sample Size | Folder |
|----------|---------------------|-------------|---------|
| 1 | Uniform(-3, 3) | n = 500 | `uniform_Z/n500/` |
| 2 | Uniform(-3, 3) | n = 1000 | `uniform_Z/n1000/` |
| 3 | N(0, 1) | n = 500 | `normal_Z/n500/` |
| 4 | N(0, 1) | n = 1000 | `normal_Z/n1000/` |

**Note**: For Normal(0,1) distribution, approximately 99.7% of Z values fall within [-3, 3], ensuring fair comparison with Uniform scenarios on the same evaluation grid.

## Folder Structure

```
Model12_SB_simulations/
├── README.md                              (this file)
├── uniform_Z/
│   ├── SB_simulation_uniform.R            (main simulation script)
│   ├── n500/
│   │   ├── outputs/
│   │   │   ├── Pointwise MSE n=500.png
│   │   │   ├── Pointwise bias n =500.png
│   │   │   └── results.json
│   │   └── results.json
│   └── n1000/
│       ├── outputs/
│       └── results.json
└── normal_Z/
    ├── SB_simulation_normal.R             (to be created)
    ├── n500/
    │   ├── outputs/
    │   └── results.json
    └── n1000/
        ├── outputs/
        └── results.json
```

## Output Files

Each scenario produces the following outputs in `nXXX/outputs/`:

- **`Pointwise MSE n=XXX.png`** - Pointwise MSE curves for all 4 functions (2×2 panel plot)
- **`Pointwise bias n =XXX.png`** - Pointwise squared bias curves for all 4 functions (2×2 panel plot)
- **`results.json`** - Numerical summary in JSON format containing:
  - `sample_size`: Sample size (n)
  - `replications`: Number of Monte Carlo replications
  - `f1`, `f2`, `f3`, `f4`: For each function:
    - `ibs`: Integrated squared bias
    - `imse`: Integrated mean squared error

## Metrics

- **IBS (Integrated Squared Bias)**: Mean of squared bias across evaluation grid
  - Formula: `mean((mean(f.hat) - f.true)²)`

- **IMSE (Integrated Mean Squared Error)**: Integrated pointwise MSE
  - Formula: `(max(x_grid) - min(x_grid)) × mean(pointwise_MSE)`

## Usage

### Running a simulation

1. Open the appropriate script (e.g., `uniform_Z/SB_simulation_uniform.R`)
2. Set working directory with `setwd()` to the desired output folder:
   - For n=500: `setwd('.../uniform_Z/n500/outputs')`
   - For n=1000: `setwd('.../uniform_Z/n1000/outputs')`
3. Set desired sample size: `n <- 500` or `n <- 1000`
4. Run the script - outputs will be saved in the working directory

### Script structure

The simulation script:
- Uses `setwd()` to determine output location
- Runs parallel Monte Carlo with `nsim = 100` replications using all available cores
- Estimates coefficient functions using wsbackfit with automatic bandwidth selection
- Computes IBS and IMSE metrics for each of the 4 functions
- Generates two 2×2 panel plots (MSE and Bias)
- Saves numerical results in `results.json`

### JSON output format

```json
{
  "sample_size": 500,
  "replications": 100,
  "f1": {
    "ibs": 0.123456,
    "imse": 0.567890
  },
  "f2": {
    "ibs": 0.234567,
    "imse": 0.678901
  },
  "f3": {
    "ibs": 0.345678,
    "imse": 0.789012
  },
  "f4": {
    "ibs": 0.456789,
    "imse": 0.890123
  }
}
```

## Expected Research Questions

1. **Sample size effect**: Does IMSE decrease as n increases (500 → 1000)?
2. **Distribution effect**: How do Uniform vs Normal Z affect estimation performance?
3. **Function-specific patterns**: Do smoother functions (e.g., f₁, f₃) show lower bias than less smooth ones (e.g., f₄)?
4. **Interaction effects**: Does the sample size effect differ between Uniform and Normal distributions?

## Date Created

2025-12-01
