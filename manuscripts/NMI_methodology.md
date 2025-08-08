## Title
Bayesian Network Meta‑Interpolation: Hierarchical, Robust, and Nonlinear Extensions for Population Transportability

## Abstract
Network Meta‑Interpolation (NMI) provides a strategy to interpolate treatment effects into target populations without assuming shared effect modifiers. We develop Bayesian NMI variants: BH‑NMI (hierarchical correlations), RB‑NMI (robustness to missingness/heavy tails), BMA‑NMI (model averaging), and GP‑NMI (Gaussian‑process nonlinearity). These methods propagate uncertainty across BLUP‑style enrichment, regression‑based interpolation, and NMA stages, improving transportability and inference reliability.

## 1. Introduction
NMI addresses settings where shared effect‑modifier assumptions are questionable. The frequentist form under‑propagates uncertainty and restricts correlation structures. Bayesian variants enable hierarchical pooling, proper uncertainty propagation, and flexible function spaces.

## 2. Methods
We formalize a three‑stage Bayesian pipeline: (1) hierarchical BLUP‑style enrichment with subgroup‑level random effects; (2) regression/interpolation to target covariate profiles with priors on slope structures; (3) NMA with generated‑quantities contrasts. GP‑NMI replaces linear interpolation with kernels, providing smooth nonlinearity.

## 3. Robustness and Sensitivity
RB‑NMI incorporates mixture components for heavy tails and missing‑data models with Bayesian imputation. Sensitivity analyses quantify residual bias under plausible violations.

## References
Harari O, et al. (2023). Network meta‑interpolation. JRSS or related venue (update final citation).

## Acknowledgments
Large language models—Claude Sonnet 4, Grok 4, ChatGPT‑3o, and ChatGPT‑5—assisted in concept refinement, drafting, and code iteration. Human authors ensured methodological rigor and validation.

