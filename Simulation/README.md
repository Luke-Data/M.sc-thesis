# Varying-Coefficient Models Simulation Study

Repository per lo studio di simulazione Monte Carlo sui modelli a coefficienti variabili (Varying-Coefficient Models) utilizzando diverse tecniche di stima e strutture di dipendenza.

## Struttura della Repository

```
Simulation/
├── README.md                           (questo file)
├── Model12_SB_simulations/             (simulazioni Smooth Backfitting e P-splines)
│   ├── README.md                       (documentazione dettagliata)
│   ├── Simulation _script/             (script R principali)
│   │   ├── Smooth_BF_VC.R             (Smoothed Backfitting)
│   │   ├── Spline_VC.R                (P-splines con GAM)
│   │   └── Permutation_test.R         (test di permutazione)
│   ├── Smooth_BF/                     (risultati Smooth Backfitting)
│   │   ├── 500_independence/          (struttura indipendente)
│   │   ├── 500_dependence/            (struttura dipendente)
│   │   └── DPG_graph.png              (grafo del modello grafico gaussiano)
│   └── P-spline/                      (risultati P-splines)
│       ├── 500_idependence/           (struttura indipendente)
│       ├── 500_dependence/            (struttura dipendente)
│       └── DPG_graph.png              (grafo del modello grafico gaussiano)
└── KLASSO simulation/                  (simulazioni Kernel LASSO)
    ├── KLASSO 3 different DGP.R       (script simulazione)
    └── Output/
        ├── Readme.md
        ├── 2_DGP_f3(Z)_non_relev.png  (risultati DGP 2)
        └── 3_DPG_f4(Z)_non_relev.png  (risultati DGP 3)
```

## Metodi di Stima Implementati

### 1. Smoothed Backfitting (wsbackfit)

**Script**: [Smooth_BF_VC.R](Model12_SB_simulations/Simulation _script/Smooth_BF_VC.R)

**Descrizione**: Algoritmo iterativo per la stima di modelli additivi generalizzati con coefficienti variabili utilizzando kernel smoothing.

**Caratteristiche**:
- Selezione automatica della bandwidth (h = -1)
- Monte Carlo con 100 replicazioni
- Parallelizzazione con `doParallel`
- Valutazione su griglia: x_grid = seq(-3, 3, length = n)

**Modelli Valutati**:

#### Struttura di Indipendenza
- **Modello**: Y = X₂f₁(X₁) + X₆f₂(X₃) + X₈f₃(X₁) + X₅f₄(X₃) + ε
- **Grafo**: p = 9 variabili, graph = "random", prob = 0.23
- **Sample size**: n = 500

#### Struttura di Dipendenza
- **Modello**: Y = X₂f₁(X₁) + X₆f₂(X₃) + X₈f₃(X₁) + X₅f₄(X₃) + ε
- **Grafo**: p = 9 variabili, graph = "random", prob = 0.4
- **Sample size**: n = 500
- **Caratteristica**: Tutte le variabili (coefficienti e smoothing) provengono dal modello grafico gaussiano

### 2. P-splines con GAM (mgcv)

**Script**: [Spline_VC.R](Model12_SB_simulations/Simulation _script/Spline_VC.R)

**Descrizione**: Stima dei modelli a coefficienti variabili utilizzando penalized splines attraverso il package `mgcv`.

**Caratteristiche**:
- Basis: P-splines (bs = 'ps')
- Metodo di smoothing: REML
- Monte Carlo con 100 replicazioni
- Stessa struttura di dipendenza del Smoothed Backfitting per confrontabilità

**Modelli Valutati**: Identici a quelli del Smoothed Backfitting per permettere confronto diretto delle performance.

### 3. Kernel LASSO (KLASSO)

**Script**: [KLASSO 3 different DGP.R](KLASSO simulation/KLASSO 3 different DGP.R)

**Descrizione**: Metodo di selezione delle variabili per modelli a coefficienti variabili tramite penalizzazione LASSO kernel-based.

**Caratteristiche**:
- Selezione automatica bandwidth tramite cross-validation (np::npscoefbw)
- Griglia di penalità: λ ∈ [0.1, 200] con 30 valori
- Algoritmo iterativo con tolleranza 1e-06
- Evaluation su 3 diversi Data Generating Processes (DGP)
