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
- Two-step estimation with different bandwidths for different functions

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
    └── Model12_SB_simulations/                  (Model 1.2 - Smooth Backfitting simulations)
        ├── README.md
        ├── Simulation _script/
        │   ├── SB_simulation.R                  (Independence structure)
        │   ├── SB_simulation_interact.R         (Dependence structure)
        │   └── SB_Graph_simulation.R            (Graph structure)
        ├── Indipendence_structure/
        │   ├── uniform_Z/
        │   │   ├── n500/
        │   │   └── n1000/
        │   └── normal_Z/
        │       ├── n500/
        │       └── n1000/
        ├── Dependece_structure/
        │   ├── n500/
        │   └── n1000/
        └── Graph_structure/
            ├── n500/
            └── 500d/
```

## Key Difference Between Models

| Aspect | Model 1.1 | Model 1.2 |
|--------|-----------|-----------|
| **Smoothing variable** | Single Z for all coefficients | Different Zⱼ for each coefficient |
| **Flexibility** | All functions vary together | Each function varies independently |
| **Estimation methods** | NW, LLR, LLP, KLASSO, Two-step | Smooth Backfitting |
| **Complexity** | Lower dimensional smoothing | Higher dimensional, requires backfitting |
