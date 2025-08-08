## Title
Bayesian Simulated Treatment Comparison: Scale‑Coherent Anchoring, Network Extensions, Robustness, and Adaptivity

## Abstract
Simulated Treatment Comparison (STC) fits an outcome model on IPD and predicts outcomes for AgD populations. Foundational critiques emphasize scale dependence and under‑propagation of uncertainty. We propose Bayesian Hierarchical STC (BH‑STC) that constructs anchored comparisons on the linear predictor scale, Network STC (N‑STC) that enforces network consistency, Robust STC (R‑STC) with measurement‑error/missing‑data models and sensitivity analysis, and Adaptive STC (A‑STC) with nonparametric components and model averaging. We formalize anchored/unanchored estimands, derive likelihoods for binary/continuous outcomes, and provide practical diagnostics.

## Methods
For binary outcomes with logit link, anchored B vs C in population P is
\[ d_{BC}(P) = [\eta_C(P) - \eta_A(P)] - [\eta_B(P) - \eta_A(P)], \]
estimated with AgD constraints on arms A and C and IPD‑based prediction of \(\eta_B(P)\) and \(\eta_A(P)\). Uncertainty from outcome model parameters and AgD sampling is fully propagated in a single Bayesian model. Network and robust extensions are described with corresponding priors and likelihoods.

## References
Ishak KJ, et al. (2015). PharmacoEconomics.\
Phillippo DM (2019). PhD Thesis, University of Bristol.

## Acknowledgments
LLMs—Claude Sonnet 4, Grok 4, ChatGPT‑3o, and ChatGPT‑5—were involved in code drafting and documentation; human authors validated and refined the methodology and results.

