# Graph Structure Simulations

## Overview

This folder contains Monte Carlo simulations evaluating the performance of the Smoothed Backfitting algorithm under **specific graphical model configurations**. The simulations assess how the dependency structure among predictors (generated via Gaussian Graphical Models) affects the estimation accuracy of varying-coefficient functions.

## Simulation Design

### Model Specification

**Statistical Model**:
```
Y = X₁·f₁(Z₁) + X₆·f₂(Z₂) + X₈·f₃(Z₁) + X₅·f₄(Z₂) + ε
```

where:
- **ε ~ N(0, 1)**: Independent error term
- **X₁, X₅, X₆, X₈**: Coefficient variables from Gaussian Graphical Model
- **Z₁, Z₂**: Smoothing variables (using X₉ and X₃ from the graphical model)

### Key Features

- **Graphical Model**: Gaussian graphical model with p = 15 predictors
  - Graph type: "random"
  - Edge probability: prob = 0.23
  - Generated using `huge.generator()` from the `huge` package
- **Shared Smoothing Variables**: Two functions share the same smoothing variable
  - f₁ and f₃ both depend on Z₁ (X₉)
  - f₂ and f₄ both depend on Z₂ (X₃)
- **Coefficient Functions**:
  - f₁(z) = z (linear)
  - f₂(z) = 2 tanh(z) (hyperbolic tangent)
  - f₃(z) = 2 exp(-0.5z²) (Gaussian-like)
  - f₄(z) = 2.5 sin(2z) exp(-0.2z²) (oscillatory with decay)

### Fixed Simulation Parameters

| Parameter | Value |
|-----------|-------|
| Sample size (n) | 500 |
| Monte Carlo replications | 100 |
| Evaluation grid | x_grid = seq(-3, 3, length = 500) |
| Bandwidth selection | Automatic (h = -1) |
| Random seed | 12 |
| Estimation method | Smoothed Backfitting (`wsbackfit::sback`) |

## Folder Structure

```
Graph_structure/
├── README.md                         (this file)
├── n500/
│   ├── EF_500.png                    (Estimated Functions plot)
│   ├── MSE_500.png                   (Pointwise MSE plot)
│   ├── bias_500.png                  (Pointwise Bias plot)
│   └── results.json                  (Numerical summary)
└── 500d/                             (Alternative configuration)
    ├── Estimate_function_500d.png    (Estimated Functions plot)
    ├── MSE_500d.png                  (Pointwise MSE plot)
    ├── Bias_500d.png                 (Pointwise Bias plot)
    └── results.json                  (Numerical summary)
```

## Two Configurations

### Configuration 1: n500/
- **Smoothing variables**: X₉ (Z₁) and X₃ (Z₂)
- **Model formula**: `Y ~ sb(X9,by=X1) + sb(X3,by=X6) + sb(X9,by=X8) + sb(X3,by=X5)`
- **Purpose**: Standard graph structure with two shared smoothing variables

### Configuration 2: 500d/
- **Smoothing variables**: X₁₁ (Z₁) and X₁₃ (Z₂)
- **Model formula**: `Y ~ sb(X11,by=X1) + sb(X13,by=X6) + sb(X11,by=X8) + sb(X13,by=X5)`
- **Purpose**: Alternative variable selection from the same graphical model
- **Note**: "d" likely indicates "different" variable indices

## Output Files

Each configuration folder contains:

### 1. Estimated Functions Plot
- **File**: `EF_*.png` or `Estimate_function_*.png`
- **Content**: 2×2 panel plot comparing true vs estimated functions
- **Features**:
  - True function (solid teal line)
  - Estimated function (dashed magenta line)
  - 95% confidence bands (gray ribbon)
  - Individual panels for each of the 4 functions

### 2. Pointwise MSE Plot
- **File**: `MSE_*.png`
- **Content**: 2×2 panel plot showing pointwise Mean Squared Error
- **Color scheme**: Blue lines with shaded area
- **Y-axis**: MSE values (range 0-10)
- **X-axis**: Evaluation grid [-3, 3]

### 3. Pointwise Bias Plot
- **File**: `bias_*.png` or `Bias_*.png`
- **Content**: 2×2 panel plot showing pointwise bias
- **Color scheme**: Green lines with shaded area
- **Y-axis**: Bias values (range approximately -0.6 to 0.6)
- **X-axis**: Evaluation grid [-3, 3]

### 4. Numerical Results (JSON)
- **File**: `results.json`
- **Format**:
```json
{
  "sample_size": 500,
  "replications": 100,
  "f1_X9X1": {
    "ibs": 0.123456,
    "imse": 0.567890
  },
  "f2_X3X6": {
    "ibs": 0.234567,
    "imse": 0.678901
  },
  "f3_X9X8": {
    "ibs": 0.345678,
    "imse": 0.789012
  },
  "f4_X3X5": {
    "ibs": 0.456789,
    "imse": 0.890123
  }
}
```
**Note**: For configuration 500d/, the keys are `f1_X11X1`, `f2_X13X6`, `f3_X11X8`, `f4_X13X5`.

## Performance Metrics

### Integrated Squared Bias (IBS)
- **Definition**: Mean of squared bias across evaluation grid
- **Formula**: `IBS = mean((E[f̂] - f_true)²)`
- **Interpretation**: Measures systematic deviation from true function
- **Lower is better**

### Integrated Mean Squared Error (IMSE)
- **Definition**: Integral of pointwise MSE over evaluation domain
- **Formula**: `IMSE = 6 × mean(pointwise_MSE)`
- **Factor 6**: Represents range width (3 - (-3) = 6)
- **Interpretation**: Overall estimation accuracy (bias² + variance)
- **Lower is better**

## Running the Simulation

### Script Location
[SB_Graph_simulation.R](../Simulation _script/SB_Graph_simulation.R)

### Steps to Run

1. **Open R/RStudio** and load the script

2. **Set working directory** (modify line 13 or 403):
   ```r
   # For n500 configuration:
   setwd('path/to/Graph_structure/n500')

   # For 500d configuration:
   setwd('path/to/Graph_structure/500d')
   ```

3. **Ensure required packages are installed**:
   ```r
   install.packages(c("foreach", "doParallel", "huge", "wsbackfit",
                      "mvtnorm", "ggplot2", "gridExtra", "patchwork",
                      "cowplot", "jsonlite", "qgraph"))
   ```

4. **Run the script**:
   - The script runs two simulations sequentially
   - First simulation: lines 1-394 (n500 configuration)
   - Second simulation: lines 396-783 (500d configuration)

5. **Output files** will be saved automatically in the working directory

### Computational Requirements

- **Parallelization**: Uses all available CPU cores minus 1
- **Memory**: Approximately 500 MB RAM per simulation
- **Runtime**: ~10-30 minutes per configuration (depends on CPU)
- **Cores used**: Detected automatically via `parallel::detectCores()`

## Research Questions

This simulation setup addresses:

1. **Graphical Model Impact**: How does the dependency structure among predictors (defined by the graphical model) affect estimation accuracy?

2. **Shared Smoothing Variables**: What is the impact of having multiple functions share the same smoothing variable (f₁ and f₃ share Z₁; f₂ and f₄ share Z₂)?

3. **Configuration Comparison**: Do different variable selections from the same graphical model (n500 vs 500d) lead to different estimation performance?

4. **Function Complexity**: How does function shape (linear, hyperbolic, Gaussian-like, oscillatory) interact with the graphical structure?

## Key Differences from Other Structures

| Aspect | Independence Structure | Dependence Structure | **Graph Structure** |
|--------|----------------------|---------------------|-------------------|
| **Predictor dependency** | Graphical model | Graphical model | **Graphical model** |
| **Smoothing variables** | Independent | From graphical model | **From graphical model** |
| **Shared smoothing** | No | Yes (complex) | **Yes (pairs)** |
| **Number of functions** | 4 varying-coefficient | 6 (4 VC + 2 smooth) | **4 varying-coefficient** |
| **Model complexity** | Moderate | High | **Moderate-High** |

## Technical Details

### Data Generation Process
1. Generate covariance matrix Σ using `huge.generator()` with random graph structure
2. Sample predictors: **X** ~ MVN(0, Σ) with p = 15 variables
3. Generate response: Y = X₁f₁(Z₁) + X₆f₂(Z₂) + X₈f₃(Z₁) + X₅f₄(Z₂) + ε
4. Estimate using: `sback(Y ~ sb(Z1,by=X1) + sb(Z2,by=X6) + sb(Z1,by=X8) + sb(Z2,by=X5))`

### Estimation Details
- **Algorithm**: Smoothed Backfitting with automatic bandwidth selection
- **Bandwidth**: h = -1 (triggers automatic selection)
- **Centering**: All estimated functions are centered: f̂(z) - mean(f̂(z))
- **Prediction**: Evaluated on a fine grid of 500 points from -3 to 3

### Error Handling
The simulation includes `tryCatch()` to handle potential estimation failures:
- Failed replications are excluded from analysis
- Final results report number of valid replications
- Matrices are resized to accommodate only valid results

## Date Created
2025-12-01

## Last Updated
2025-12-05
