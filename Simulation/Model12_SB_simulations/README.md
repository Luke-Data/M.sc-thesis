# Model 1.2 Varying-Coefficient Simulations

## Overview

Studio di simulazione Monte Carlo per confrontare due metodi di stima per modelli a coefficienti variabili sotto diverse strutture di dipendenza generate da Gaussian Graphical Models (GGM).

### Metodi Comparati

1. **Smoothed Backfitting (wsbackfit)**: Algoritmo iterativo kernel-based con selezione automatica della bandwidth
2. **P-splines (GAM/mgcv)**: Penalized splines con smoothing parameter selection via REML

### Parametri Comuni delle Simulazioni

- **Monte Carlo replications**: nsim = 100
- **Sample size**: n = 500
- **Evaluation grid**: x_grid = seq(-3, 3, length = n)
- **Data generation**: Multivariate normal distributions con struttura GGM (via `huge.generator`)
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

**Interpretazione**: Se p-value è basso, rifiutiamo H₀ → la funzione varia effettivamente con z₂.

### Setup del Test

```r
n = 500
nsim = 100
f1(x) = 2·exp(-0.5x²)
f2(x) = 2.5·sin(2x)·exp(-0.2x²)

x1 ~ Beta(2, 3)
x2 ~ Poisson(5)
z1, z2 ~ N(0, 1)
ε ~ N(0, 1)

Y = x₁·f₁(z₁) + x₂·f₂(z₂) + ε
```

## Esecuzione delle Simulazioni

### Prerequisiti R

```r
install.packages(c(
  "wsbackfit",      # Smooth backfitting
  "mgcv",           # GAM e P-splines
  "huge",           # Gaussian Graphical Models
  "mvtnorm",        # Multivariate normal
  "foreach",        # Parallel loops
  "doParallel",     # Parallel backend
  "ggplot2",        # Visualizzazione
  "patchwork",      # Layout multipli
  "cowplot",        # Grafici combinati
  "qgraph"          # Visualizzazione grafi
))
```

### Procedura per Smooth Backfitting

1. Aprire [Smooth_BF_VC.R](Simulation _script/Smooth_BF_VC.R)

2. Modificare il `setwd()` alla riga 12:
   ```r
   # Per indipendenza:
   setwd('path/to/Model12_SB_simulations/Smooth_BF/500_independence')

   # Per dipendenza:
   setwd('path/to/Model12_SB_simulations/Smooth_BF/500_dependence')
   ```

3. Verificare/modificare i parametri GGM:
   ```r
   n <- 500
   p <- 9  # numero variabili

   # Per indipendenza: prob = 0.23
   Gr <- huge.generator(n, p, graph="random", prob=0.23)

   # Per dipendenza: prob = 0.4
   Gr <- huge.generator(n, p, graph="random", prob=0.4)
   ```

4. Eseguire lo script completo

5. Output salvati automaticamente in working directory

### Procedura per P-splines

1. Aprire [Spline_VC.R](Simulation _script/Spline_VC.R)

2. Modificare il `setwd()` alla riga 11:
   ```r
   # Per indipendenza:
   setwd('path/to/Model12_SB_simulations/P-spline/500_idependence')

   # Per dipendenza:
   setwd('path/to/Model12_SB_simulations/P-spline/500_dependence')
   ```

3. Verificare i parametri (identici a Smooth BF):
   ```r
   n <- 500
   p <- 9
   Gr <- huge.generator(n, p, graph="random", prob=0.23)  # o 0.4
   ```

4. Eseguire lo script completo

### Procedura per Test di Permutazione

1. Aprire [Permutation_test.R](Simulation _script/Permutation_test.R)

2. Lo script contiene due sezioni separate:
   - Righe 1-71: Test 1 (permutazione coefficiente)
   - Righe 72-107: Test 2 (permutazione smoothing variable)

3. Eseguire la sezione desiderata o entrambe

4. I p-value vengono stampati in console

5. Grafici generati interattivamente

## Caratteristiche Tecniche degli Script

### Smoothed Backfitting (Smooth_BF_VC.R)

**Metodo di stima**:
```r
fit <- sback(Y ~ sb(X1, by=X2, h=-1) + sb(X3, by=X6, h=-1) +
                 sb(X1, by=X8, h=-1) + sb(X3, by=X5, h=-1),
             data=data)
```

**Caratteristiche**:
- `h = -1`: bandwidth selection automatica
- `sb()`: smooth backfitting term
- `by=`: specifica la variabile coefficiente
- Predizione su griglia con `predict.sback()`

**Estrazione funzioni stimate**:
```r
f.hat1 <- pred$coeff['X2'] +
          pred$coeff['X2:X1'] * x_grid +
          pred$peffects[,1]
f.hat1 <- f.hat1 - mean(f.hat1)  # centratura
```

### P-splines (Spline_VC.R)

**Metodo di stima**:
```r
fit <- gam(Y ~ s(X1, bs='ps', by=X2) + s(X3, bs='ps', by=X6) +
               s(X1, bs='ps', by=X8) + s(X3, bs='ps', by=X5),
           data=data, method="REML")
```

**Caratteristiche**:
- `bs='ps'`: P-splines basis
- `method="REML"`: Restricted Maximum Likelihood per smoothing parameter
- `by=`: variabile coefficiente (interazione)
- Predizione: `predict(fit, type="terms")`

**Estrazione funzioni stimate**:
```r
f.hat1 <- pred[, grep("s\\(X1\\):X2", colnames(pred))]
f.hat1 <- f.hat1 - mean(f.hat1)  # centratura
```

### Parallelizzazione

Entrambi gli script usano:
```r
n_cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

sim.results <- foreach(i = 1:nsim,
                       .packages = c(...),
                       .export = c(...)) %dopar% {
  # Codice simulazione
}

stopCluster(cl)
```

## Domande di Ricerca

1. **Confronto Metodi**: Smoothed Backfitting vs P-splines
   - Quale metodo ha IBS e IMSE più bassi?
   - La performance relativa dipende dalla forma funzionale?
   - Quale metodo è più robusto alla dipendenza tra variabili?

2. **Effetto Struttura di Dipendenza**:
   - Come cambia la performance quando prob passa da 0.23 a 0.4?
   - La dipendenza tra variabili aiuta o danneggia la stima?
   - Quale metodo è meno sensibile alla struttura di dipendenza?

3. **Difficoltà delle Funzioni**:
   - Quale tra f₁ (lineare), f₂ (sigmoide), f₃ (gaussiana), f₄ (oscillante) è più difficile da stimare?
   - L'ordinamento di difficoltà è lo stesso per entrambi i metodi?
   - Il bias è più alto per funzioni oscillanti?

4. **Bandwidth vs Smoothing Parameter**:
   - La selezione automatica della bandwidth (wsbackfit) è comparabile a REML (mgcv)?
   - Ci sono differenze sistematiche nel grado di smoothing applicato?

## Note Tecniche

### Riproducibilità
- Seed globale per GGM: `set.seed(12)`
- Seed per ogni replicazione MC: `set.seed(i)` dentro il loop `%dopar%`
- Questo garantisce risultati identici su diverse macchine

### Gestione Errori
Entrambi gli script hanno `tryCatch` per gestire errori di stima:
```r
tryCatch({
  # codice stima
}, error = function(e) {
  return(NULL)  # replicazione fallita
})
```

### Centratura delle Funzioni
Tutte le funzioni (vere e stimate) sono centrate:
```r
f.true <- f(x_grid) - mean(f(x_grid))
f.hat <- f.hat - mean(f.hat)
```

Questo elimina l'intercetta e permette confronto diretto delle forme funzionali.

### Griglia di Valutazione
```r
x_grid <- seq(-3, 3, length = n)  # n = 500 punti
```
- Range [-3, 3] copre la maggior parte della massa di probabilità per variabili ~ N(0,1)
- 500 punti garantiscono risoluzione fine per curve smooth
- Stessa griglia per vero e stimato → confronto diretto

## Risultati Preliminari

Basati sui JSON disponibili:

### Smoothed Backfitting
**Indipendenza** (prob=0.23):
- f₁: IBS=0.001, IMSE=3.476
- f₂: IBS=0.006, IMSE=1.125
- f₃: IBS=0.015, IMSE=5.637
- f₄: IBS=0.117, IMSE=2.458

**Dipendenza** (prob=0.4):
- f₁: IBS=0.002, IMSE=0.808
- f₂: IBS=0.010, IMSE=1.152
- f₃: IBS=0.011, IMSE=0.953
- f₄: IBS=0.145, IMSE=1.756

**Osservazione**: La dipendenza riduce IMSE per f₁ e f₃, ma aumenta leggermente per f₂ e f₄.

### P-splines (GAM)
**Indipendenza**:
- f₁: IBS=0.000, IMSE=0.062
- f₂: IBS=0.001, IMSE=0.315
- f₃: IBS=0.007, IMSE=0.333
- f₄: IBS=0.005, IMSE=1.011

**Osservazione**: P-splines mostra IMSE drammaticamente più bassi, specialmente per f₁, f₂, f₃. Suggerisce maggiore efficienza in questo setup.

## Date

- **Creazione repository**: 2024-12-01
- **Ultima esecuzione simulazioni**: 2024-12-11
- **Ultimo aggiornamento documentazione**: 2025-12-11

## Riferimenti

- Mammen, E., Linton, O., & Nielsen, J. (1999). The existence and asymptotic properties of a backfitting projection algorithm under weak conditions. *Annals of Statistics*, 27(5), 1443-1490.
- Wood, S. N. (2017). *Generalized Additive Models: An Introduction with R* (2nd ed.). CRC Press.
- Hastie, T., & Tibshirani, R. (1993). Varying-coefficient models. *Journal of the Royal Statistical Society: Series B*, 55(4), 757-779.
