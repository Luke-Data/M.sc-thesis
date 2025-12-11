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

**Data Generating Processes**:

#### DGP 1 (Baseline)
```
Y = X₁·sin(Z) + X₂·4Z(1-Z) + ε
```
- n = 400
- X₁ ~ N(2, 1.5²), X₂ ~ Beta(0.5, 0.5), X₃ ~ Gamma(0.5, 1)
- Z ~ Uniform(0, 1)
- Tutte le funzioni sono rilevanti

#### DGP 2 (f₃ non rilevante)
```
Y = X₁·f₁(Z) + X₂·f₂(Z) + X₃·0 + ε
```
- f₃(Z) = cos(πZ) viene annullata
- **Obiettivo**: Verificare capacità di identificare funzioni nulle

#### DGP 3 (f₄ non rilevante)
```
Y = X₁·f₁(Z) + X₂·f₂(Z) + X₃·f₃(Z) + X₄·0 + ε
```
- f₄(Z) non rilevante
- **Obiettivo**: Test di robustezza con variabili spurie

## Funzioni Coefficiente Utilizzate

Le simulazioni utilizzano un set comune di funzioni test:

```r
f1(z) = z                                    # Lineare
f2(z) = 2·tanh(z)                           # Sigmoide
f3(z) = 2·exp(-0.5z²)                       # Gaussiana
f4(z) = 2.5·sin(2z)·exp(-0.2z²)            # Oscillante con decadimento
```

Per KLASSO:
```r
f1(z) = sin(z)                              # Sinusoidale
f2(z) = 4z(1-z)                             # Quadratica
f3(z) = cos(πz)                             # Cosinusoidale
```

## Strutture di Dipendenza

### Modello Grafico Gaussiano

Tutte le simulazioni utilizzano Gaussian Graphical Models (GGM) generati tramite `huge.generator`:

- **Tipo di grafo**: Random
- **Parametri di connessione**:
  - Indipendenza: prob = 0.23 (15 variabili)
  - Dipendenza: prob = 0.4 (9 variabili)

### Struttura di Indipendenza
Le variabili coefficiente (Xᵢ) sono generate dal GGM, ma le variabili di smoothing rimangono indipendenti tra le funzioni.

### Struttura di Dipendenza
Sia le variabili coefficiente che le variabili di smoothing condividono la struttura di dipendenza del GGM.

## Metriche di Performance

Tutte le simulazioni calcolano:

### IBS (Integrated Squared Bias)
```r
IBS = mean((mean(f.hat) - f.true)²)
```
Misura il bias quadratico medio integrato sulla griglia di valutazione.

### IMSE (Integrated Mean Squared Error)
```r
IMSE = 6 × mean(pointwise_MSE)
```
Dove il fattore 6 rappresenta l'ampiezza dell'intervallo: max(x_grid) - min(x_grid) = 3 - (-3) = 6.

## Grafici Generati

### Smooth Backfitting e P-splines

Ogni configurazione produce:

1. **Pointwise MSE**: MSE punto per punto per tutte le funzioni (pannelli 2×2)
   - File: `Pointwise_MSE_*.png`

2. **Pointwise Bias**: Bias punto per punto per tutte le funzioni (pannelli 2×2)
   - File: `Pointwise_bias_*.png`

3. **Estimated Functions**: Confronto funzioni vere vs stimate
   - File: `VC_estimate_*.png`

4. **DGP Graph**: Visualizzazione del grafo del modello gaussiano
   - File: `DPG_graph.png`

**Esempi**:
- [Pointwise_MSE_SB.png](Model12_SB_simulations/Smooth_BF/500_independence/Pointwise_MSE_SB.png)
- [VC_estimate_P-spline_d.png](Model12_SB_simulations/P-spline/500_dependence/VC_estimate_P-spline_d.png)

### KLASSO

1. **Variable Selection Results**: Path dei coefficienti al variare di λ
   - [2_DGP_f3(Z)_non_relev.png](KLASSO simulation/Output/2_DGP_f3(Z)_non_relev.png)
   - [3_DPG_f4(Z)_non_relev.png](KLASSO simulation/Output/3_DPG_f4(Z)_non_relev.png)

## Test di Permutazione

**Script**: [Permutation_test.R](Model12_SB_simulations/Simulation _script/Permutation_test.R)

Implementa due tipi di test di permutazione per verificare la rilevanza delle funzioni coefficiente:

### Test 1: Permutazione della variabile coefficiente
- **H₀**: f₂(z₂) = 0 (fluttuazioni attorno allo zero)
- Permuta X₂ mantenendo Z₂ fisso
- Confronta MSE del modello originale vs MSE dei modelli permutati
- nsim = 100 permutazioni

### Test 2: Permutazione della variabile smoothing
- **H₀**: f(z₂) = costante (non varia con z₂)
- Permuta Z₂ mantenendo X₂ fisso
- Confronta MSE del modello originale vs MSE dei modelli permutati
- nsim = 100 permutazioni

**Setup**:
```r
n = 500
f1(x) = 2·exp(-0.5x²)
f2(x) = 2.5·sin(2x)·exp(-0.2x²)
Y = x₁·f₁(z₁) + x₂·f₂(z₂) + ε
```

**Output**: p-value empirico basato sulla proporzione di MSE permutati ≤ MSE originale.

## File JSON di Output

Ogni configurazione di simulazione salva i risultati numerici in formato JSON:

```json
{
  "sample_size": 500,
  "replications": 100,
  "f1_X2X1": {
    "ibs": 0.0123456,
    "imse": 0.0567890
  },
  "f2_X6X3": {
    "ibs": 0.0234567,
    "imse": 0.0678901
  },
  "f3_X8X1": {
    "ibs": 0.0345678,
    "imse": 0.0789012
  },
  "f4_X5X3": {
    "ibs": 0.0456789,
    "imse": 0.0890123
  }
}
```

**Naming convention**:
- `results.json` (Smooth Backfitting)
- `results_gam.json` (P-splines)
- `results_sb.json` (Smooth Backfitting dipendenza)

## Esecuzione delle Simulazioni

### Requisiti

```r
# Pacchetti richiesti
install.packages(c(
  "wsbackfit",      # Smooth backfitting
  "mgcv",           # GAM e P-splines
  "huge",           # Generazione GGM
  "mvtnorm",        # Distribuzione normale multivariata
  "np",             # Nonparametric kernel regression
  "foreach",        # Parallel loops
  "doParallel",     # Backend parallelizzazione
  "ggplot2",        # Visualizzazione
  "patchwork",      # Composizione grafici
  "qgraph"          # Visualizzazione grafi
))
```

### Procedura Generale

1. **Aprire lo script desiderato**
2. **Modificare il working directory** per l'output desiderato
3. **Impostare il sample size** (se applicabile): `n <- 500` o `n <- 1000`
4. **Eseguire lo script completo**
5. **I risultati verranno salvati** nella working directory impostata

### Esempi Specifici

#### Smooth Backfitting - Indipendenza
```r
setwd('path/to/Model12_SB_simulations/Smooth_BF/500_independence')
n <- 500
# Esegui Smooth_BF_VC.R
```

#### P-splines - Dipendenza
```r
setwd('path/to/Model12_SB_simulations/P-spline/500_dependence')
n <- 500
# Esegui Spline_VC.R
```

#### KLASSO
```r
# Modifica i DGP commentando/decommentando le sezioni appropriate
# Esegui KLASSO 3 different DGP.R
```

## Domande di Ricerca

1. **Confronto Metodi**: Come si confrontano Smooth Backfitting e P-splines in termini di IBS e IMSE?

2. **Effetto Dipendenza**: Qual è l'impatto della struttura di dipendenza sulla qualità della stima?

3. **Variable Selection**: Il KLASSO identifica correttamente le funzioni non rilevanti?

4. **Robustezza**: Come variano le performance al variare della forma funzionale (lineare, gaussiana, oscillante)?

5. **Bandwidth Selection**: L'effetto della selezione automatica della bandwidth confrontata con approcci fissi.

## Note Tecniche

- **Parallelizzazione**: Tutti gli script Monte Carlo utilizzano `doParallel` con `n_cores - 1` core
- **Seed**: I seed sono impostati per riproducibilità (seed = 12 per GGM, seed = i per ogni replicazione MC)
- **Centratura**: Tutte le funzioni stimate sono centrate: `f.hat <- f.hat - mean(f.hat)`
- **Griglia di valutazione**: Uniforme su [-3, 3] con `length = n` punti

## Data e Versioni

- **Creazione**: 2024-08-25 (KLASSO)
- **Ultimo aggiornamento**: 2025-12-11 (Smooth BF & P-splines)
- **R version**: ≥ 4.0.0 consigliata

## Riferimenti

Per dettagli specifici sui singoli esperimenti, consultare:
- [Model12_SB_simulations/README.md](Model12_SB_simulations/README.md) - Documentazione dettagliata Smooth Backfitting e P-splines
- [KLASSO simulation/Output/Readme.md](KLASSO simulation/Output/Readme.md) - Risultati KLASSO
