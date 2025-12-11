# M.sc Thesis - Varying Coefficient Models

This repository contains the implementation and simulations for varying coefficient models research.

## Models

### Model 1.1
**All coefficient functions share the same smoothing variable**

E(Y|X=x,Z=z) = x₁ · f₁(z) + x₂ · f₂(z) + ... + xₐ · fₐ(z)

where all coefficient functions f₁(z), f₂(z), ..., fₐ(z) are functions of the **same smoothing variable Z**.

**Location**: [Model11/](Model11/)

**Methods implemented**:
- NW (Nadaraya-Watson): VCM with constant local fit
- LLR (Local Linear Regression): VCM with local linear function fit
- LLP (Local Polynomial): VCM with local cubic fit
- KLASSO: Group lasso penalty for nonparametric regression
- MSE bandwidth selector for LLR
- two-step estimation with different bandwidths for different functions

### Model 1.2
**Each coefficient function can have its own smoothing variable**

E(Y|X=x,Z=z) = x₁ · f₁(z₁) + x₂ · f₂(z₂) + ... + xₐ · fₐ(zₐ)

where each coefficient function fⱼ(zⱼ) can depend on a **different smoothing variable Zⱼ**.

**Location**: [Model12/](Model12/)

**Method implemented**:
- Smooth Backfitting algorithm

## Repository Structure

```
M.sc-thesis/
├── README.md                                    (this file)
│
├── Model11/                                     (Model 1.1: same smoothing variable)
│   ├── README.md
│   ├── KLASSO.R
│   ├── Model 1.1.R
│   ├── mse_vcm.R
│   ├── two-step.R
│   └── Plot/
│       ├── README.md
│       ├── custom KLASSO trial/
│       └── Custom LLR & mgcv package estimate/
│
├── Model12/                                     (Model 1.2: different smoothing variables)
│   ├── README.md
│   └── Model 1.2 - backfitting.R
│
└── Simulation/
    ├── KLASSO simulation/                       (Model 1.1 - KLASSO simulations)
    │   └── Output/
    │
    ├── LLR_MSE_2step simulation/                (Model 1.1 - LLR and Two-step simulations)
    │   └── Output/
    │
    └── Model12_SB_simulations/                  (Model 1.2 - Varying-Coefficient comparison)
        ├── README.md
        ├── Simulation _script/
        │   ├── Smooth_BF_VC.R                   (Smoothed Backfitting method)
        │   ├── Spline_VC.R                      (P-splines/GAM method)
        │   └── Permutation_test.R               (Permutation tests)
        ├── Smooth_BF/                           (Smoothed Backfitting results)
        │   ├── DPG_graph.png
        │   ├── 500_independence/
        │   │   ├── Pointwise_MSE_SB.png
        │   │   ├── Pointwise_bias_SB.png
        │   │   ├── VC_estimate_SB.png
        │   │   └── results.json
        │   └── 500_dependence/
        │       ├── Pointwise_MSE_SB_d.png
        │       ├── Pointwise_bias_SB_d.png
        │       ├── VC_estimate_SB_d.png
        │       └── results_sb.json
        └── P-spline/                            (P-splines results)
            ├── DPG_graph.png
            ├── 500_idependence/
            │   ├── Pointwise_MSE_P-spline.png
            │   ├── Pointwise_bias_P-spline.png
            │   ├── VC estimate P-spline.png
            │   └── results_gam.json
            └── 500_dependence/
                ├── Pointwise_MSE_P-spline_d.png
                ├── Pointwise_bias_P-spline_d.png
                ├── VC_estimate_P-spline_d.png
                └── results_gam.json
```

## Key Difference Between Models

| Aspect | Model 1.1 | Model 1.2 |
|--------|-----------|-----------|
| **Smoothing variable** | Single Z for all coefficients | Different Zⱼ for each coefficient |
| **Flexibility** | All functions vary together | Each function varies independently |
| **Estimation methods** | NW, LLR, LLP, KLASSO, Two-step | Smooth Backfitting, P-splines (GAM) |
| **Complexity** | Lower dimensional smoothing | Higher dimensional, requires backfitting |

## Model 1.2 Simulations - Method Comparison

The Model 1.2 simulations compare two estimation approaches for varying-coefficient models:

### Estimation Methods

1. **Smoothed Backfitting (wsbackfit)**
   - Iterative kernel-based algorithm
   - Automatic bandwidth selection (h = -1)
   - Implemented via `wsbackfit` R package

2. **P-splines (GAM/mgcv)**
   - Penalized splines with basis functions
   - Smoothing parameter selection via REML
   - Implemented via `mgcv::gam()` with `s()` smoothers

### Dependency Structures Tested

#### Independence Structure
- **Purpose**: Smoothing variables (Z) are independent from other model variables
- **Model**: Y = X₂·f₁(Z₁) + X₆·f₂(Z₂) + X₈·f₃(Z₁) + X₅·f₄(Z₂) + ε
- **Folders**: `Smooth_BF/500_independence/`, `P-spline/500_idependence/`

#### Dependence Structure
- **Purpose**: Both coefficient and smoothing variables share dependency structure (Gaussian Graphical Model)
- **Model**: Y = X₂·f₁(X₁) + X₆·f₂(X₃) + X₈·f₃(X₁) + X₅·f₄(X₃) + ε
- **Folders**: `Smooth_BF/500_dependence/`, `P-spline/500_dependence/`

### Common Simulation Parameters

- **Sample size**: n = 500
- **Monte Carlo replications**: 100
- **Data generation**: Multivariate normal with Gaussian Graphical Model (p = 9, prob = 0.4)
- **True coefficient functions**:
  - f₁(z) = z (linear)
  - f₂(z) = 2·tanh(z) (hyperbolic tangent)
  - f₃(z) = 2·exp(-0.5z²) (Gaussian-like)
  - f₄(z) = 2.5·sin(2z)·exp(-0.2z²) (oscillatory with decay)

### Performance Metrics

- **IBS (Integrated Squared Bias)**: `mean((mean(f.hat) - f.true)²)`
- **IMSE (Integrated Mean Squared Error)**: `6 × mean(pointwise_MSE)`

### Permutation Tests

The repository also includes permutation test implementations for testing:
1. **Coefficient variable relevance**: H₀: f(z) = 0
2. **Smoothing variable effect**: H₀: f(z) = constant

See [Simulation/Model12_SB_simulations/README.md](Simulation/Model12_SB_simulations/README.md) for detailed documentation.
