
# **epi-econ-causal**

**Replication and simulation code for the working paper:**
**“A Unified Causal Inference Framework Bridging Epidemiology and Econometrics”**

**Authors:**
Awa Diop, Ababacar S. Gueye, Cheikh I. Nokho, Kossi D. Abalo,
Janick Weberpals, Nima Hejazi, Denis Talbot

---

## Overview

This repository contains simulation and estimation code supporting a methodological study that bridges causal inference approaches used in **epidemiology** and **econometrics**.

The goal of the project is to compare identification strategies and estimators across disciplines within a unified framework, highlighting both shared foundations and key conceptual differences.

The repository includes:

* Monte Carlo simulation studies
* Data-generating processes (DGPs) tailored to different identification designs
* Implementations of standard causal estimators
* Reproducible scripts for generating tables and figures

---

## Simulation Designs

Two primary simulation settings are included:

###  **Unmeasured Confounding Scenario**

File:
`sim_unmeasured_confounding.R`

* The DGP includes an **unobserved confounder (U)** that affects both treatment and outcome
* (U) is **not included** in estimation models
* Evaluates robustness of estimators under unmeasured confounding

**Important:**
The treatment effect is constant across individuals. As a result, the **LATE** estimated by the IV design is numerically close to the **ATT**. This occurs due to the simulation design and does **not** hold in general.

---

### **No Unmeasured Confounding Scenario**

File:
`sim_no_unmeasured_confounding.R`

* All confounders are observed and included in adjustment
* Serves as a benchmark where identification assumptions hold

---

## 📊 Estimators Compared

| Identification Strategy   | Estimator                    |
| ------------------------- | ---------------------------- |
| Confounding control (ATT) | Matching                     |
|                           | IPTW                         |
|                           | G-computation                |
|                           | Doubly Robust (AIPTW)        |
| Panel methods (ATT)       | Difference-in-Differences    |
| Instrumental Variables    | IV (LATE)                    |
| Regression Discontinuity  | RDD (local effect at cutoff) |

**Note:**
DiD, IV, and RDD are evaluated under **separate data-generating designs** aligned with their identifying assumptions, rather than being applied to the same cross-sectional ATT DGP.

---

## Running the Simulations

Run either script directly in R:

```r
source("sim_unmeasured_confounding.R")
source("sim_no_unmeasured_confounding.R")
```

Each script:

* Runs Monte Carlo simulations
* Produces summary performance metrics (Bias, RMSE, Coverage)
* Generates bias distribution plots

---

## 📌 Purpose

This repository is designed as a **transparent and reproducible companion** to the methodological paper. It illustrates how different causal inference traditions approach identification and estimation under varying assumptions.

---

## 📜 License

MIT
