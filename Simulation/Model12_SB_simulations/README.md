# Model 1.2 Smoothed Backfitting Simulations

## Overview

Monte Carlo study evaluating wsbackfit performance under three different dependence structures using Gaussian Graphical Models.

### Common Fixed Elements

- **Monte Carlo replications**: nsim = 100
- **Estimation method**: wsbackfit with automatic bandwidth selection (h = -1)
- **Evaluation grid**: x_grid = seq(-3, 3, length = n)
- **Sample sizes**: n = 500, n = 1000
- **Data generation**: Multivariate normal distributions with graphical model structure (using `huge.generator`)

## Three Simulation Structures

### 1. Independence Structure

**Purpose**: Evaluate performance when coefficient variables (Xᵢ) are generated from a Gaussian Graphical Model, but smoothing variables remain independent across functions.

**Script**: `SB_simulation.R`

**Model**: Y = X₁f₁(Z₁) + X₆f₂(Z₂) + X₈f₃(Z₁) + X₅f₄(Z₂) + ε

- **Graphical Model**: p = 15 predictors, graph = "random", prob = 0.23
- **True coefficient functions**:
  - f₁(z) = z
  - f₂(z) = 2 tanh(z)
  - f₃(z) = 2 exp(-0.5z²)
  - f₄(z) = 2.5 sin(2z) exp(-0.2z²)

**Folder structure**:
```
Indipendence_structure/
├── uniform_Z/
│   ├── n500/      (Z ~ Uniform(-3, 3))
│   └── n1000/
└── normal_Z/
    ├── n500/      (Z ~ N(0, 1))
    └── n1000/
```

### 2. Dependence Structure (Interaction)

**Purpose**: Evaluate performance when both coefficient variables and smoothing variables share dependencies, with additional non-varying smooth terms.

**Script**: `SB_simulation_interact.R`

**Model**: Y = X₁f₁(X₂) + X₆f₂(X₃) + X₁₀f₃(X₈) + X₅f₄(X₃) + f₅(X₂) + f₃(X₇) + ε

- **Graphical Model**: p = 10 predictors, graph = "random", prob = 0.4
- **True coefficient functions**:
  - f₁(z) = z
  - f₂(z) = 2 tanh(z)
  - f₃(z) = 2 exp(-0.5z²)
  - f₄(z) = 0.1z³
  - f₅(z) = 2.5 sin(2z) exp(-0.2z²)
- **Note**: Includes 6 functions total (4 varying-coefficient + 2 smooth terms)

**Folder structure**:
```
Dependece_structure/
├── n500/
└── n1000/
```

### 3. Graph Structure

**Purpose**: Evaluate performance under a specific graphical model configuration with both dependent and independent smoothing variables.

**Script**: `SB_Graph_simulation.R`

**Model**: Y = X₁f₁(Z₁) + X₆f₂(Z₂) + X₈f₃(Z₃) + X₅f₄(Z₄) + ε

- **Graphical Model**: p = 15 predictors, graph = "random", prob = 0.23
- **True coefficient functions**: Same as Independence Structure

**Folder structure**:
```
Graph_structure/
├── n500/
└── 500d/     (alternative configuration)
```

## Folder Structure

```
Model12_SB_simulations/
├── README.md                              (this file)
├── Simulation _script/
│   ├── SB_simulation.R                    (Independence structure)
│   ├── SB_simulation_interact.R           (Dependence structure with interactions)
│   └── SB_Graph_simulation.R              (Graph structure)
├── Indipendence_structure/
│   ├── uniform_Z/
│   │   ├── n500/
│   │   │   ├── Pointwise MSE n=500.png
│   │   │   ├── Pointwise bias n =500.png
│   │   │   └── results.json
│   │   └── n1000/
│   │       ├── Poinwise MSE n=1000.png
│   │       ├── Poinwise Bias n=1000.png
│   │       └── results.json
│   └── normal_Z/
│       ├── n500/                          (empty - to be generated)
│       └── n1000/
│           ├── Poinwise MSE normal n=1000.png
│           ├── Pointwise bias normal n=1000.png
│           ├── Norm_est_func_n=1000.png
│           └── results.json
├── Dependece_structure/
│   ├── n500/
│   │   ├── EF_interaction_n=500.png
│   │   ├── MSE interaction n=500.png
│   │   ├── Bias interaction n=500.png
│   │   └── results.json
│   └── n1000/
│       ├── MSE_interaction_n1000.png
│       ├── Bias_interaction_n1000.png
│       ├── true_f_interaction_n1000.png
│       └── results.json
└── Graph_structure/
    ├── n500/
    │   ├── EF_500.png
    │   ├── MSE_500.png
    │   ├── bias_500.png
    │   └── results.json
    └── 500d/
        ├── Estimate_function_500d.png
        ├── MSE_500d.png
        ├── Bias_500d.png
        └── results.json
```

## Output Files

### Independence Structure
Each scenario folder contains:
- **Pointwise MSE plot** - MSE curves for all 4 functions (2×2 panel plot)
- **Pointwise Bias plot** - Bias curves for all 4 functions (2×2 panel plot)
- **Estimated Functions plot** - True vs estimated functions comparison
- **results.json** - Numerical summary with IBS and IMSE for each function

### Dependence Structure
Each sample size folder contains:
- **MSE interaction plot** - MSE curves for all 6 functions (2×3 panel plot)
- **Bias interaction plot** - Bias curves for all 6 functions (2×3 panel plot)
- **EF/True function plot** - True vs estimated functions comparison
- **results.json** - Numerical summary with IBS and IMSE for all 6 functions

### Graph Structure
Each folder contains:
- **MSE plot** - MSE curves for all 4 functions
- **Bias plot** - Bias curves for all 4 functions
- **Estimated Functions plot** - True vs estimated functions
- **results.json** - Numerical summary with IBS and IMSE for each function

## Metrics

- **IBS (Integrated Squared Bias)**: Mean of squared bias across evaluation grid
  - Integral approximation formula: `mean((mean(f.hat) - f.true)²)`

- **IMSE (Integrated Mean Squared Error)**: Integrated pointwise MSE
  - Integral approximation formula: `6 × mean(pointwise_MSE)`
  - Factor 6 represents the range width: max(x_grid) - min(x_grid) = 3 - (-3) = 6

## Common Script Features

All simulation scripts:
- Run parallel Monte Carlo with `nsim = 100` replications using all available cores (`doParallel`)
- Estimate coefficient functions using `wsbackfit` with automatic bandwidth selection (h = -1)
- Compute IBS and IMSE metrics for each function
- Generate visualization plots using `ggplot2` and `patchwork`
- Save numerical results to JSON format

## JSON Output Format

### Independence and Graph Structure
```json
{
  "sample_size": 500,
  "replications": 100,
  "f1_XiXj": {
    "ibs": 0.123456,
    "imse": 0.567890
  },
  "f2_XkXl": {
    "ibs": 0.234567,
    "imse": 0.678901
  },
  "f3_XmXn": {
    "ibs": 0.345678,
    "imse": 0.789012
  },
  "f4_XpXq": {
    "ibs": 0.456789,
    "imse": 0.890123
  }
}
```

### Dependence Structure (6 functions)
```json
{
  "sample_size": 500,
  "replications": 100,
  "f1_X1X2": { "ibs": ..., "imse": ... },
  "f2_X6X3": { "ibs": ..., "imse": ... },
  "f3_X10X8": { "ibs": ..., "imse": ... },
  "f4_X5X3": { "ibs": ..., "imse": ... },
  "f5_X2": { "ibs": ..., "imse": ... },
  "f3_X7": { "ibs": ..., "imse": ... }
}
```

## Date Created

2025-12-01

## Last Updated

2025-12-05
