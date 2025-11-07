# JMNCC: Joint Modeling of Longitudinal Biomarker and Survival Outcomes with Competing Risk in Nested Case-Control Studies
## Description
**JMNCC** is an R implementation of a joint modeling framework for **longitudinal biomarkers** and **survival outcomes** in the presenece of **competing risks** under a **nested case-control (NCC) design**.

This package provides two complementary estimation strategies:
1. **fJM-NCC** - Full likelihood-based inference using all observed longitudinal data and survival outcome, treating the unobserved longitudinal measurements as missing at random.
2. **wJM-NCC** - Inverse probability weighting likelihood-based inference using only data from NCC sub-cohort, accounting for the potential selection bias in the NCC sampling.

## Model Framework

Let $T_i$ denote the survival time, $\delta_i \in (0, 1, \dots, K)$ the event indicator for censoring or competing event $k$, and $Y_i(t)$ the longitudinal biomarker for subject $i$.  Using an indicator variable $R_i=1$ for subject $i$ included in the NCC sub-cohort, $R_i=0$ otherwise.

The JM-NCC model consists of two linked sub-models:

1. **Longitudinal sub-model**'

   A generalized linear mixed-effects model capturing the temporal trajectory of biomarker measurements:

   $$g(E(Y_i (t)|b_i))=\gamma^T X^{(1)}_i (t)  + b^T_i Z_i ,$$

   where $b_i \sim N(0, D)$ are random effects.
   
2. **Survival sub-model (cause-specific hazards)**

   For competing event $k$:

   $$\lambda_{k}(t|b_i)=\lambda_{0k} (t) \exp(\beta_k \cdot g^{-1}( \gamma^T X^{(1)}_i (t)  + b^T_i Z_i ) + \alpha_k^T \cdot X^{(2)}_i),$$

   where $\lambda_{0k} (t)$ is approximated by a piecewise-constant baseline hazard, and $g^{-1}( \gamma^T X^{(1)}_i (t)+ b^T_i Z_i)$ is the subject-specific fitted trajectory from the longitudinal model.

## Likelihood Formulation

### 1.fJM-NCC (Full likelihood)

The fJM-NCC method incorporates all observed data and survival outcome, treating the unobserved longitudinal outcomes as missing at random:

$$l_{full} = \sum^N_{i=1; R_i=0} log \int f(T_i, \delta_i|b_i, X_i) f(b_i) d b_i+ \sum^N_{i=1; R_i=1} log \int f(T_i, \delta_i|b_i, X_i) f(Y_i|b_i, X_i) f(b_i) d b_i.$$
 
This approach is statistically efficient but computationally intensive, as it integrates over both the longitudinal and survival components for all subjects.

### 2.wJM-NCC (Inverse probability weighting likelihood)

The wJM-NCC method constructs an inverse probability weighting likelihood function using only data from the NCC sub-cohort:
 
$$l_{wt} = \sum^N_{i=1; R_i=1} w_i log \int f(T_i, \delta_i|b_i, X_i) f(Y_i|b_i, X_i) f(b_i) d b_i,$$

where $w_i=(1-\prod_{l \in S_i} (1-\frac{m_l}{n_l}) )^{-1}$ presents the inverse probability of inclusion as a control.

### Component Definitions

- $f(T_i, \delta_i|b_i, X_i)$: cause-specific model with piecewise-constant baseline hazard and association parameter $\beta_k$,
- $f(Y_i|b_i, X_i)$: longitudinal model,
- $f(b_i)$: random effects distribution (multivariate normal)

Integration over $b_i$ is performed using **Gauss–Hermite quadrature**, and the survival part is integrated with **Gauss–Kronrod quadrature** for time-dependent effects.

## Implementation and Features

- Supports **competing risks** 
- Compatible with 'optim()' for maximum likelihood estimation
- Provides **gradient functions** for efficient optimization
- Handles **missing longitudinal data** in fJM-NCC method

## Installation
```r
# Install development version from GitHub
# (requires devtools)
install.packages("devtools")
devtools::install_github("Zhaoyn-oss/JMNCC")
```

# JMNCC
