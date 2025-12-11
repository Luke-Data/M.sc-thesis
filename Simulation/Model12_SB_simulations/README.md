# Model 1.2 Varying-Coefficient Simulations

## Overview

Studio di simulazione Monte Carlo per confrontare due metodi di stima per modelli a coefficienti variabili sotto diverse strutture di dipendenza 

### Metodi Comparati

1. **Smoothed Backfitting (wsbackfit)**: Algoritmo iterativo kernel-based con selezione automatica della bandwidth
2. **P-splines (GAM/mgcv)**: Penalized splines con smoothing parameter selection via REML

### Parametri Comuni delle Simulazioni

- **Monte Carlo replications**: nsim = 100
- **Sample size**: n = 500
- **Evaluation grid**: x_grid = seq(-3, 3, length = n)
- **Data generation**: Multivariate normal distributions 
- **Parallelization**: `doParallel` con tutti i core disponibili - 1

## Strutture di Dipendenza Analizzate

### 1. Struttura di Indipendenza

**Obiettivo**: Valutare le performance quando le variabili coefficiente (Xᵢ) sono generate dal GGM, ma le variabili di smoothing sono indipendenti tra le funzioni.

**Script**:
- Smoothed Backfitting: [Smooth_BF_VC.R](Simulation _script/Smooth_BF_VC.R)
- P-splines: [Spline_VC.R](Simulation _script/Spline_VC.R)

**Modello**:
```
Y = X₂·f₁(X₁) + X₆·f₂(X₃) + X₈·f₃(X₁) + X₅·f₄(X₃) + ε
```

**Data Generating Process**:
- **Graphical Model**: p = 9 variabili, graph = "random", prob = 0.23
- **Funzioni coefficiente**:
  - f₁(z) = z
  - f₂(z) = 2·tanh(z)
  - f₃(z) = 2·exp(-0.5z²)
  - f₄(z) = 2.5·sin(2z)·exp(-0.2z²)

**Note**: In questo setup, X₁ e X₃ (variabili di smoothing) possono condividere dipendenze attraverso il GGM, ma sono usate in funzioni diverse. Le variabili coefficiente (X₂, X₆, X₈, X₅) provengono anch'esse dal GGM.

**Output folders**:
```
Smooth_BF/500_independence/
P-spline/500_idependence/
```

### 2. Struttura di Dipendenza

**Obiettivo**: Valutare le performance quando sia le variabili coefficiente che le variabili di smoothing condividono la struttura di dipendenza del GGM.

**Script**:
- Smoothed Backfitting: [Smooth_BF_VC.R](Simulation _script/Smooth_BF_VC.R)
- P-splines: [Spline_VC.R](Simulation _script/Spline_VC.R)

**Modello**:
```
Y = X₂·f₁(X₁) + X₆·f₂(X₃) + X₈·f₃(X₁) + X₅·f₄(X₃) + ε
```

**Data Generating Process**:
- **Graphical Model**: p = 9 variabili, graph = "random", prob = 0.4 (maggiore connettività)
- **Funzioni coefficiente**: Identiche alla struttura di indipendenza
- **Caratteristica distintiva**: Maggiore probabilità di connessioni nel grafo (0.4 vs 0.23)

**Output folders**:
```
Smooth_BF/500_dependence/
P-spline/500_dependence/
```

### Differenza Chiave tra le Strutture

La differenza principale tra le due strutture non è nel modello, ma nella **probabilità di connessione del GGM**:
- **Indipendenza**: prob = 0.23 → grafo più sparso → minore dipendenza tra variabili
- **Dipendenza**: prob = 0.4 → grafo più denso → maggiore dipendenza tra variabili

Questo permette di studiare come l'intensità delle dipendenze nel GGM influenzi le performance di stima.

## Struttura delle Cartelle

```
Model12_SB_simulations/
├── README.md                           (questo file)
├── Simulation _script/
│   ├── Smooth_BF_VC.R                 (Smoothed Backfitting)
│   ├── Spline_VC.R                    (P-splines con GAM)
│   └── Permutation_test.R             (Test di permutazione)
├── Smooth_BF/                         (Risultati Smoothed Backfitting)
│   ├── DPG_graph.png                  (Visualizzazione GGM)
│   ├── results.json                   (Risultati vecchi - deprecato)
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
└── P-spline/                          (Risultati P-splines)
    ├── DPG_graph.png                  (Visualizzazione GGM)
    ├── 500_idependence/               (nota: typo nel nome folder)
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

## Output Generati

Ogni configurazione (metodo × struttura di dipendenza) produce:

### 1. Pointwise MSE Plot
File: `Pointwise_MSE_*.png`

Pannelli 2×2 che mostrano il Mean Squared Error punto per punto per ciascuna delle 4 funzioni coefficiente lungo la griglia di valutazione [-3, 3].

**Interpretazione**:
- Curve più basse = migliore performance
- Permette di identificare regioni dello spazio dove la stima è più difficile
- Confronto visivo tra le 4 diverse forme funzionali

### 2. Pointwise Bias Plot
File: `Pointwise_bias_*.png`

Pannelli 2×2 che mostrano il bias punto per punto (media delle stime - funzione vera) per ciascuna funzione.

**Interpretazione**:
- Oscillazioni attorno allo zero = bias basso
- Deviazioni sistematiche indicano distorsione della stima
- Identifica se il metodo tende a sovrastimare o sottostimare in certe regioni

### 3. Estimated Functions Plot
File: `VC_estimate_*.png` o `VC estimate *.png`

Pannelli 2×2 che sovrappongono:
- Funzione vera (linea colorata/nera)
- Funzioni stimate da ogni replicazione MC (linee grigie trasparenti)
- Media delle funzioni stimate (linea blu/rossa spessa)

**Interpretazione**:
- Permette di visualizzare la variabilità delle stime
- La fascia di linee grigie rappresenta la distribuzione delle stime
- Confronto immediato tra funzione vera e stima media

### 4. DGP Graph
File: `DPG_graph.png` (uno per metodo)

Visualizzazione del Gaussian Graphical Model utilizzato per generare i dati.

**Caratteristiche**:
- Nodi = variabili (X₁, X₂, ..., X₉)
- Archi = dipendenze condizionali
- Generato con `qgraph`
- Colori: edge.color="#7A6F8E", color="#EDE8F2"

### 5. Results JSON
File: `results.json`, `results_sb.json`, `results_gam.json`

Sommario numerico delle performance con IBS e IMSE per ogni funzione.

## Metriche di Performance

### IBS (Integrated Squared Bias)

**Formula**:
```r
IBS = mean((mean(f.hat) - f.true)²)
```

**Interpretazione**:
- Misura il bias quadratico medio integrato
- Quantifica la distorsione sistematica della stima
- Valori più bassi = minore bias
- Non dipende dalla variabilità delle stime

### IMSE (Integrated Mean Squared Error)

**Formula**:
```r
IMSE = 6 × mean(pointwise_MSE)
```

**Interpretazione**:
- Misura l'errore quadratico medio integrato
- Combina bias e varianza
- Il fattore 6 = range width: max(x_grid) - min(x_grid) = 3 - (-3)
- Metrica principale per confrontare metodi
- Valori più bassi = migliore performance complessiva

### Esempio JSON Output

#### Smoothed Backfitting - Indipendenza
```json
{
  "sample_size": 500,
  "replications": 100,
  "f1_X1Z1": {
    "ibs": 0.001,
    "imse": 3.476
  },
  "f2_X6Z2": {
    "ibs": 0.006,
    "imse": 1.125
  },
  "f3_X8Z1": {
    "ibs": 0.015,
    "imse": 5.637
  },
  "f4_X5Z2": {
    "ibs": 0.117,
    "imse": 2.458
  }
}
```

#### P-splines - Dipendenza
```json
{
  "sample_size": 500,
  "replications": 100,
  "method": "GAM P-splines (REML)",
  "f1_X2X1": {
    "ibs": 0.002,
    "imse": 0.808
  },
  "f2_X6X3": {
    "ibs": 0.010,
    "imse": 1.152
  },
  "f3_X8X1": {
    "ibs": 0.011,
    "imse": 0.953
  },
  "f4_X5X3": {
    "ibs": 0.145,
    "imse": 1.756
  }
}
```

## Test di Permutazione

**Script**: [Permutation_test.R](Simulation _script/Permutation_test.R)

Implementa due approcci per testare la rilevanza delle funzioni coefficiente usando Smoothed Backfitting.

### Test 1: Permutazione della Variabile Coefficiente

**H₀**: f₂(z₂) = 0 (la funzione oscilla attorno allo zero senza effetto sistematico)

**Procedura**:
1. Stima il modello completo: `Y ~ sb(z₁, by=x₁) + sb(z₂, by=x₂)`
2. Calcola MSE del modello: `mse.l`
3. Per i = 1 a nsim:
   - Permuta casualmente x₂ → x₂.p
   - Stima `Y ~ sb(z₁, by=x₁) + sb(z₂, by=x₂.p)`
   - Calcola MSE permutato: `mse.b[i]`
4. p-value = proporzione di `mse.b ≤ mse.l`

**Interpretazione**: Se p-value è basso, rifiutiamo H₀ → la funzione è significativa.

### Test 2: Permutazione della Variabile di Smoothing

**H₀**: f(z₂) = costante (la funzione non varia con z₂)

**Procedura**:
1. Stima il modello completo: `Y ~ sb(z₁, by=x₁) + sb(z₂, by=x₂)`
2. Calcola MSE del modello: `mse.l`
3. Per i = 1 a nsim:
   - Permuta casualmente z₂ → z₂.p
   - Stima `Y ~ sb(z₁, by=x₁) + sb(z₂.p, by=x₂)`
   - Calcola MSE permutato: `mse.b[i]`
4. p-value = proporzione di `mse.b ≤ mse.l`

