## Title
Bayesian Component Network Meta‑Analysis with Interactions, Adaptive Selection, and Hierarchical Heterogeneity

## Abstract
Component network meta‑analysis (cNMA) decomposes treatments into components to estimate component‑level effects and their combinations. We propose three Bayesian advances—BCI‑cNMA (interaction modeling), ACS‑cNMA (adaptive component selection), and HCE‑cNMA (hierarchical component‑specific heterogeneity)—to address non‑additivity, model uncertainty, and heterogeneity propagation. The methods are implemented in Stan with cmdstanr, include WAIC/LOO diagnostics, and are validated via simulation under additive, synergistic, antagonistic, and sparse networks. We provide formal models, priors, inference strategies, and guidance for application.

## 1. Introduction
Component network meta‑analysis enables modular inference when interventions comprise multiple elements. Existing approaches often rely on additivity assumptions that may be violated in practice and inadequately model component‑level heterogeneity and model uncertainty. We present a Bayesian framework that (i) flexibly captures interactions, (ii) selects predictive components with uncertainty, and (iii) models component‑specific heterogeneity hierarchically.

## 2. Model Formulation
Let treatment \( k \) consist of components in set \( C_k \subseteq \{1,\ldots, J\} \). For study \( s \) and arm \( a \in k \), define linear predictor
\[ \eta_{sa} = \mu_s + \sum_{j\in C_k} \alpha_j + \sum_{(j,\ell)\in C_k\times C_k, j<\ell} \alpha_{j\ell} + \epsilon_{sa}, \]
with baseline \( \mu_s \sim \mathcal{N}(0, \sigma_\mu^2) \). The additive terms \( \alpha_j \) are component main effects; the pair terms \( \alpha_{j\ell} \) represent interactions. Outcome‑specific likelihoods follow GLM families, e.g., Bernoulli‑logit or Normal.

### 2.1 BCI‑cNMA (Interactions)
Assign hierarchical priors to interaction terms,
\[ \alpha_{j\ell} \sim \mathcal{N}(0, \tau_{int}^2), \quad \tau_{int} \sim \mathrm{half\text{-}Normal}(0, \sigma_{int}), \]
enabling shrinkage towards additivity while retaining power to detect synergy/antagonism.

### 2.2 ACS‑cNMA (Adaptive Selection)
Use spike‑and‑slab or horseshoe priors for main/interaction effects, e.g.,
\[ \alpha_j \sim \pi_j \mathcal{N}(0, \sigma_1^2) + (1-\pi_j)\, \delta_0, \quad \pi_j \sim \mathrm{Beta}(a,b), \]
with analogous structure for \( \alpha_{j\ell} \), to select predictive components and interactions.

### 2.3 HCE‑cNMA (Hierarchical Heterogeneity)
Model component‑specific heterogeneity across studies: for study‑specific deviations \( \delta_{sj} \),
\[ \alpha_{j}^{(s)} = \alpha_j + \delta_{sj}, \quad \delta_{sj} \sim \mathcal{N}(0, \tau_{j}^2), \quad \tau_j \sim \mathrm{half\text{-}Normal}(0, \sigma_\tau), \]
allowing components to vary in their between‑study variance.

## 3. Inference and Computation
Models are fit with HMC/NUTS via Stan. We compute posterior contrasts for any pair of treatment combinations and perform model comparison using WAIC/LOO. Posterior predictive checks assess fit to study‑arm outcomes.

## 4. Practical Guidance
We recommend ACS‑cNMA when the component set is large or sparsely used, BCI‑cNMA when synergy/antagonism is plausible, and HCE‑cNMA when between‑study variation plausibly differs by component. Priors should reflect expected sparsity and plausible effect magnitudes.

## References
Beliveau A, et al. (2017). Component network meta‑analysis. Research Synthesis Methods.\
Rücker G, Petropoulou M, Schwarzer G (2019). Network meta‑analysis of multicomponent interventions. Statistics in Medicine.\
Efthimiou O, et al. (2022). Methods for cNMA with multiple outcomes. Research Synthesis Methods.

## Acknowledgments
We acknowledge the use of large language models as co‑development tools in ideation, drafting, and code generation: Claude Sonnet 4, Grok 4, ChatGPT‑3o, and ChatGPT‑5. Human authors curated assumptions, validated models, designed simulations, and ensured scientific correctness.

