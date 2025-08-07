# Novel Component Network Meta-Analysis Methodologies
# Author: Research Collaboration
# Date: 2025
# Purpose: Develop and test new Bayesian cNMA methods

# Load required libraries
library(cmdstanr)
library(posterior)
library(bayesplot)
library(loo)
library(brms)
library(MCMCpack)
library(mvtnorm)
library(netmeta)
library(gemtc)
library(BUGSnet)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(knitr)
library(xtable)
library(parallel)
library(foreach)
library(doParallel)

# Set up parallel processing
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Methodology 1: Bayesian Component Interaction Network Meta-Analysis (BCI-NMA)
# =============================================================================

#' Bayesian Component Interaction Network Meta-Analysis
#' 
#' This function implements a novel Bayesian approach that models synergistic
#' and antagonistic interactions between treatment components, moving beyond
#' simple additive assumptions.
#' 
#' @param data Data frame with columns: study, treatment, n, events (for binary) or mean, sd (for continuous)
#' @param component_matrix Matrix indicating which components are in each treatment
#' @param interaction_terms Logical matrix indicating which component pairs may interact
#' @param outcome_type "binary" or "continuous"
#' @param prior_interaction_sd Standard deviation for interaction effect priors
#' @param n_chains Number of MCMC chains
#' @param n_iter Number of iterations per chain
#' @param n_warmup Number of warmup iterations
#' 
#' @return List containing posterior samples, model diagnostics, and summaries

bci_nma <- function(data, component_matrix, interaction_terms = NULL, 
                    outcome_type = "binary", prior_interaction_sd = 0.5,
                    n_chains = 4, n_iter = 4000, n_warmup = 2000) {
  
  # Prepare data for Stan
  n_studies <- length(unique(data$study))
  n_treatments <- nrow(component_matrix)
  n_components <- ncol(component_matrix)
  
  # Create study-level data
  studies <- unique(data$study)
  study_data <- data.frame(
    study_id = 1:n_studies,
    study_name = studies
  )
  
  # Merge with original data
  data <- merge(data, study_data, by.x = "study", by.y = "study_name")
  
  # Create treatment mapping
  treatments <- rownames(component_matrix)
  treatment_map <- data.frame(
    treatment_id = 1:n_treatments,
    treatment_name = treatments
  )
  
  data <- merge(data, treatment_map, by.x = "treatment", by.y = "treatment_name")
  
  # Prepare interaction matrix
  if (is.null(interaction_terms)) {
    # Default: all pairwise interactions possible
    interaction_terms <- matrix(1, n_components, n_components)
    diag(interaction_terms) <- 0
  }
  
  n_interactions <- sum(upper.tri(interaction_terms) & interaction_terms == 1)
  
  # Create interaction index
  interaction_idx <- which(upper.tri(interaction_terms) & interaction_terms == 1, arr.ind = TRUE)
  
  if (outcome_type == "binary") {
    stan_data <- list(
      N = nrow(data),
      n_studies = n_studies,
      n_treatments = n_treatments,
      n_components = n_components,
      n_interactions = n_interactions,
      study = data$study_id,
      treatment = data$treatment_id,
      events = data$events,
      n_patients = data$n,
      component_matrix = component_matrix,
      interaction_idx = interaction_idx,
      prior_interaction_sd = prior_interaction_sd
    )
    
    stan_model_code <- "
    data {
      int<lower=0> N;
      int<lower=0> n_studies;
      int<lower=0> n_treatments;
      int<lower=0> n_components;
      int<lower=0> n_interactions;
      int<lower=1, upper=n_studies> study[N];
      int<lower=1, upper=n_treatments> treatment[N];
      int<lower=0> events[N];
      int<lower=0> n_patients[N];
      matrix[n_treatments, n_components] component_matrix;
      int interaction_idx[n_interactions, 2];
      real<lower=0> prior_interaction_sd;
    }
    
    parameters {
      vector[n_studies] mu;  // Study baselines
      vector[n_components] beta;  // Component main effects
      vector[n_interactions] gamma;  // Component interaction effects
      real<lower=0> tau;  // Between-study heterogeneity
      matrix[n_studies, n_components] delta;  // Study-specific component effects
    }
    
    transformed parameters {
      vector[N] logit_p;
      
      for (i in 1:N) {
        logit_p[i] = mu[study[i]];
        
        // Add main component effects
        for (c in 1:n_components) {
          logit_p[i] += component_matrix[treatment[i], c] * 
                       (beta[c] + delta[study[i], c]);
        }
        
        // Add interaction effects
        for (int_i in 1:n_interactions) {
          int c1 = interaction_idx[int_i, 1];
          int c2 = interaction_idx[int_i, 2];
          logit_p[i] += component_matrix[treatment[i], c1] * 
                       component_matrix[treatment[i], c2] * gamma[int_i];
        }
      }
    }
    
    model {
      // Priors
      mu ~ normal(0, 2);
      beta ~ normal(0, 1);
      gamma ~ normal(0, prior_interaction_sd);
      tau ~ normal(0, 1);
      
      for (c in 1:n_components) {
        delta[c] ~ normal(0, tau);
      }
      
      // Likelihood
      events ~ binomial_logit(n_patients, logit_p);
    }
    
    generated quantities {
      vector[n_treatments] treatment_effects;
      matrix[n_treatments, n_treatments] treatment_comparisons;
      vector[N] log_lik;
      
      // Calculate treatment effects
      for (t in 1:n_treatments) {
        treatment_effects[t] = 0;
        for (c in 1:n_components) {
          treatment_effects[t] += component_matrix[t, c] * beta[c];
        }
        
        // Add interaction effects
        for (int_i in 1:n_interactions) {
          int c1 = interaction_idx[int_i, 1];
          int c2 = interaction_idx[int_i, 2];
          treatment_effects[t] += component_matrix[t, c1] * 
                                 component_matrix[t, c2] * gamma[int_i];
        }
      }
      
      // Calculate pairwise comparisons
      for (i in 1:n_treatments) {
        for (j in 1:n_treatments) {
          treatment_comparisons[i, j] = treatment_effects[i] - treatment_effects[j];
        }
      }
      
      // Calculate log likelihood for WAIC
      for (i in 1:N) {
        log_lik[i] = binomial_logit_lpmf(events[i] | n_patients[i], logit_p[i]);
      }
    }
    "
  } else {
    # Continuous outcome model
    stan_data <- list(
      N = nrow(data),
      n_studies = n_studies,
      n_treatments = n_treatments,
      n_components = n_components,
      n_interactions = n_interactions,
      study = data$study_id,
      treatment = data$treatment_id,
      y = data$mean,
      se = data$sd / sqrt(data$n),
      component_matrix = component_matrix,
      interaction_idx = interaction_idx,
      prior_interaction_sd = prior_interaction_sd
    )
    
    stan_model_code <- "
    data {
      int<lower=0> N;
      int<lower=0> n_studies;
      int<lower=0> n_treatments;
      int<lower=0> n_components;
      int<lower=0> n_interactions;
      int<lower=1, upper=n_studies> study[N];
      int<lower=1, upper=n_treatments> treatment[N];
      vector[N] y;
      vector<lower=0>[N] se;
      matrix[n_treatments, n_components] component_matrix;
      int interaction_idx[n_interactions, 2];
      real<lower=0> prior_interaction_sd;
    }
    
    parameters {
      vector[n_studies] mu;
      vector[n_components] beta;
      vector[n_interactions] gamma;
      real<lower=0> tau;
      vector[n_studies] delta[n_components];
    }
    
    transformed parameters {
      vector[N] theta;
      
      for (i in 1:N) {
        theta[i] = mu[study[i]];
        
        for (c in 1:n_components) {
          theta[i] += component_matrix[treatment[i], c] * 
                     (beta[c] + delta[study[i], c]);
        }
        
        for (int_i in 1:n_interactions) {
          int c1 = interaction_idx[int_i, 1];
          int c2 = interaction_idx[int_i, 2];
          theta[i] += component_matrix[treatment[i], c1] * 
                     component_matrix[treatment[i], c2] * gamma[int_i];
        }
      }
    }
    
    model {
      mu ~ normal(0, 2);
      beta ~ normal(0, 1);
      gamma ~ normal(0, prior_interaction_sd);
      tau ~ normal(0, 1);
      
      for (c in 1:n_components) {
        delta[c] ~ normal(0, tau);
      }
      
      y ~ normal(theta, se);
    }
    
    generated quantities {
      vector[n_treatments] treatment_effects;
      matrix[n_treatments, n_treatments] treatment_comparisons;
      vector[N] log_lik;
      
      for (t in 1:n_treatments) {
        treatment_effects[t] = 0;
        for (c in 1:n_components) {
          treatment_effects[t] += component_matrix[t, c] * beta[c];
        }
        
        for (int_i in 1:n_interactions) {
          int c1 = interaction_idx[int_i, 1];
          int c2 = interaction_idx[int_i, 2];
          treatment_effects[t] += component_matrix[t, c1] * 
                                 component_matrix[t, c2] * gamma[int_i];
        }
      }
      
      for (i in 1:n_treatments) {
        for (j in 1:n_treatments) {
          treatment_comparisons[i, j] = treatment_effects[i] - treatment_effects[j];
        }
      }
      
      // Calculate log likelihood for WAIC
      for (i in 1:N) {
        log_lik[i] = normal_lpdf(y[i] | theta[i], se[i]);
      }
    }
    "
  }
  
  # Write Stan model to temporary file
  temp_stan_file <- tempfile(fileext = ".stan")
  writeLines(stan_model_code, temp_stan_file)
  
  # Compile and fit model
  stan_model <- cmdstan_model(temp_stan_file)
  
  fit <- stan_model$sample(
    data = stan_data,
    chains = n_chains,
    iter_sampling = n_iter - n_warmup,
    iter_warmup = n_warmup,
    parallel_chains = min(n_chains, detectCores()),
    refresh = 500
  )
  
  # Extract results
  posterior_draws <- fit$draws()
  posterior <- as_draws_df(posterior_draws)
  
  # Calculate summaries using posterior package
  component_summary <- posterior::summarise_draws(posterior_draws, "beta")
  component_effects <- data.frame(
    component = paste0("Component_", 1:n_components),
    mean = component_summary$mean,
    sd = component_summary$sd,
    q2.5 = component_summary$q5,
    q97.5 = component_summary$q95
  )
  
  interaction_summary <- posterior::summarise_draws(posterior_draws, "gamma")
  interaction_effects <- data.frame(
    interaction = paste0("Interaction_", 1:n_interactions),
    mean = interaction_summary$mean,
    sd = interaction_summary$sd,
    q2.5 = interaction_summary$q5,
    q97.5 = interaction_summary$q95
  )
  
  treatment_summary <- posterior::summarise_draws(posterior_draws, "treatment_effects")
  treatment_effects <- data.frame(
    treatment = treatments,
    mean = treatment_summary$mean,
    sd = treatment_summary$sd,
    q2.5 = treatment_summary$q5,
    q97.5 = treatment_summary$q95
  )
  
  return(list(
    fit = fit,
    posterior = posterior,
    component_effects = component_effects,
    interaction_effects = interaction_effects,
    treatment_effects = treatment_effects,
    model_diagnostics = fit$summary(),
    data = stan_data
  ))
}

# Methodology 2: Adaptive Component Selection Network Meta-Analysis (ACS-NMA)
# ============================================================================

#' Adaptive Component Selection Network Meta-Analysis
#' 
#' This function uses Bayesian model selection to determine optimal component
#' decomposition strategies for complex interventions.
#' 
#' @param data Data frame with study data
#' @param possible_components List of possible component matrices to consider
#' @param outcome_type "binary" or "continuous"
#' @param prior_model_probs Prior probabilities for each component decomposition
#' 
#' @return List containing model comparison results and best model

acs_nma <- function(data, possible_components, outcome_type = "binary",
                    prior_model_probs = NULL) {
  
  n_models <- length(possible_components)
  
  if (is.null(prior_model_probs)) {
    prior_model_probs <- rep(1/n_models, n_models)
  }
  
  # Fit each possible model
  model_fits <- vector("list", n_models)
  log_likelihoods <- numeric(n_models)
  waic_values <- numeric(n_models)
  
  for (i in 1:n_models) {
    cat("Fitting model", i, "of", n_models, "\n")
    
    # Fit BCI-NMA with current component matrix
    model_fits[[i]] <- bci_nma(data, possible_components[[i]], 
                               outcome_type = outcome_type,
                               n_chains = 2, n_iter = 2000, n_warmup = 1000)
    
    # Calculate model fit statistics
    log_lik <- model_fits[[i]]$fit$draws("log_lik")
    waic_values[i] <- loo::waic(log_lik)$estimates["waic", "Estimate"]
  }
  
  # Calculate model weights using WAIC
  delta_waic <- waic_values - min(waic_values)
  waic_weights <- exp(-0.5 * delta_waic) / sum(exp(-0.5 * delta_waic))
  
  # Bayesian model averaging
  best_model_idx <- which.min(waic_values)
  
  results <- list(
    model_fits = model_fits,
    waic_values = waic_values,
    waic_weights = waic_weights,
    best_model = model_fits[[best_model_idx]],
    best_model_idx = best_model_idx,
    model_comparison = data.frame(
      model = 1:n_models,
      waic = waic_values,
      delta_waic = delta_waic,
      weight = waic_weights
    )
  )
  
  return(results)
}

# Methodology 3: Hierarchical Component Effects Network Meta-Analysis (HCE-NMA)
# ==============================================================================

#' Hierarchical Component Effects Network Meta-Analysis
#' 
#' This function accounts for component-specific heterogeneity sources using
#' hierarchical Bayesian models.
#' 
#' @param data Data frame with study data including study-level covariates
#' @param component_matrix Matrix indicating components in each treatment
#' @param study_covariates Matrix of study-level covariates
#' @param outcome_type "binary" or "continuous"
#' 
#' @return Fitted hierarchical model

hce_nma <- function(data, component_matrix, study_covariates = NULL,
                    outcome_type = "binary") {
  
  # Prepare data
  n_studies <- length(unique(data$study))
  n_treatments <- nrow(component_matrix)
  n_components <- ncol(component_matrix)
  
  if (is.null(study_covariates)) {
    study_covariates <- matrix(1, n_studies, 1)  # Intercept only
  }
  n_covariates <- ncol(study_covariates)
  
  # Prepare Stan data
  studies <- unique(data$study)
  study_data <- data.frame(study_id = 1:n_studies, study_name = studies)
  data <- merge(data, study_data, by.x = "study", by.y = "study_name")
  
  treatments <- rownames(component_matrix)
  treatment_map <- data.frame(treatment_id = 1:n_treatments, treatment_name = treatments)
  data <- merge(data, treatment_map, by.x = "treatment", by.y = "treatment_name")
  
  if (outcome_type == "binary") {
    stan_data <- list(
      N = nrow(data),
      n_studies = n_studies,
      n_treatments = n_treatments,
      n_components = n_components,
      n_covariates = n_covariates,
      study = data$study_id,
      treatment = data$treatment_id,
      events = data$events,
      n_patients = data$n,
      component_matrix = component_matrix,
      study_covariates = study_covariates
    )
    
    stan_model_code <- "
    data {
      int<lower=0> N;
      int<lower=0> n_studies;
      int<lower=0> n_treatments;
      int<lower=0> n_components;
      int<lower=0> n_covariates;
      int<lower=1, upper=n_studies> study[N];
      int<lower=1, upper=n_treatments> treatment[N];
      int<lower=0> events[N];
      int<lower=0> n_patients[N];
      matrix[n_treatments, n_components] component_matrix;
      matrix[n_studies, n_covariates] study_covariates;
    }
    
    parameters {
      vector[n_studies] mu;
      matrix[n_components, n_covariates] alpha;  // Component-covariate interactions
      vector[n_components] beta_0;  // Component main effects
      vector<lower=0>[n_components] tau_comp;  // Component-specific heterogeneity
      matrix[n_studies, n_components] delta;  // Study-component interactions
    }
    
    transformed parameters {
      vector[N] logit_p;
      matrix[n_studies, n_components] beta_study;
      
      // Calculate study-specific component effects
      for (s in 1:n_studies) {
        for (c in 1:n_components) {
          beta_study[s, c] = beta_0[c];
          for (k in 1:n_covariates) {
            beta_study[s, c] += alpha[c, k] * study_covariates[s, k];
          }
          beta_study[s, c] += delta[s, c];
        }
      }
      
      for (i in 1:N) {
        logit_p[i] = mu[study[i]];
        for (c in 1:n_components) {
          logit_p[i] += component_matrix[treatment[i], c] * beta_study[study[i], c];
        }
      }
    }
    
    model {
      // Priors
      mu ~ normal(0, 2);
      beta_0 ~ normal(0, 1);
      tau_comp ~ half_normal(0, 1);
      
      for (c in 1:n_components) {
        alpha[c] ~ normal(0, 0.5);
        delta[, c] ~ normal(0, tau_comp[c]);
      }
      
      // Likelihood
      events ~ binomial_logit(n_patients, logit_p);
    }
    
    generated quantities {
      vector[n_treatments] treatment_effects_avg;
      
      // Average treatment effects across studies
      for (t in 1:n_treatments) {
        treatment_effects_avg[t] = 0;
        for (c in 1:n_components) {
          treatment_effects_avg[t] += component_matrix[t, c] * beta_0[c];
        }
      }
    }
    "
  }
  
  # Write Stan model to temporary file
  temp_stan_file <- tempfile(fileext = ".stan")
  writeLines(stan_model_code, temp_stan_file)
  
  # Compile and fit model
  stan_model <- cmdstan_model(temp_stan_file)
  
  fit <- stan_model$sample(
    data = stan_data,
    chains = 4,
    iter_sampling = 2000,
    iter_warmup = 2000,
    parallel_chains = min(4, detectCores()),
    refresh = 500
  )
  
  return(list(
    fit = fit,
    posterior = as_draws_df(fit$draws()),
    data = stan_data
  ))
}

# Traditional/Gold Standard Methods for Comparison
# ===============================================

#' Traditional Additive Component Network Meta-Analysis
#' 
#' Implements the standard additive cNMA approach for comparison
#' 
#' @param data Study data
#' @param component_matrix Component decomposition matrix
#' @param outcome_type Type of outcome
#' 
#' @return Fitted traditional model

traditional_cnma <- function(data, component_matrix, outcome_type = "binary") {
  
  # Use netcomb from netmeta package for traditional approach
  if (outcome_type == "binary") {
    # Convert to pairwise format for netmeta
    pairwise_data <- pairwise(
      treat = data$treatment,
      event = data$events,
      n = data$n,
      studlab = data$study,
      data = data
    )
    
    # Fit standard NMA
    nma_result <- netmeta(
      TE = pairwise_data$TE,
      seTE = pairwise_data$seTE,
      treat1 = pairwise_data$treat1,
      treat2 = pairwise_data$treat2,
      studlab = pairwise_data$studlab,
      comb.random = TRUE,
      comb.fixed = FALSE
    )
    
    # Fit component model
    comp_result <- netcomb(
      nma_result,
      C.matrix = component_matrix
    )
    
    return(comp_result)
  }
}

# Utility Functions
# ================

#' Generate simulation data for cNMA
#' 
#' @param n_studies Number of studies
#' @param treatments Vector of treatment names
#' @param component_matrix Matrix defining treatment components
#' @param component_effects True component effects
#' @param interaction_effects True interaction effects
#' @param outcome_type "binary" or "continuous"
#' @param study_sizes Vector of study sizes
#' 
#' @return Simulated dataset

generate_cnma_data <- function(n_studies, treatments, component_matrix,
                               component_effects, interaction_effects = NULL,
                               outcome_type = "binary", study_sizes = NULL) {
  
  if (is.null(study_sizes)) {
    study_sizes <- sample(50:200, n_studies, replace = TRUE)
  }
  
  n_treatments <- length(treatments)
  n_components <- length(component_effects)
  
  # Generate study designs (which treatments are compared in each study)
  study_designs <- list()
  for (i in 1:n_studies) {
    n_arms <- sample(2:min(4, n_treatments), 1, prob = c(0.6, 0.3, 0.1)[1:(min(4, n_treatments)-1)])
    study_designs[[i]] <- sample(treatments, n_arms)
  }
  
  # Generate data
  data_list <- list()
  
  for (i in 1:n_studies) {
    study_treatments <- study_designs[[i]]
    n_arms <- length(study_treatments)
    
    # Calculate true treatment effects
    true_effects <- numeric(n_arms)
    for (j in 1:n_arms) {
      treatment_idx <- which(treatments == study_treatments[j])
      
      # Main component effects
      true_effects[j] <- sum(component_matrix[treatment_idx, ] * component_effects)
      
      # Add interaction effects if specified
      if (!is.null(interaction_effects)) {
        components_present <- which(component_matrix[treatment_idx, ] == 1)
        if (length(components_present) >= 2) {
          for (k1 in 1:(length(components_present)-1)) {
            for (k2 in (k1+1):length(components_present)) {
              interaction_idx <- paste0(components_present[k1], "_", components_present[k2])
              if (interaction_idx %in% names(interaction_effects)) {
                true_effects[j] <- true_effects[j] + interaction_effects[[interaction_idx]]
              }
            }
          }
        }
      }
    }
    
    # Generate study baseline
    study_baseline <- rnorm(1, 0, 0.5)
    
    # Generate arm-level data
    for (j in 1:n_arms) {
      arm_size <- round(study_sizes[i] / n_arms)
      
      if (outcome_type == "binary") {
        logit_p <- study_baseline + true_effects[j] + rnorm(1, 0, 0.2)  # Add some heterogeneity
        p <- plogis(logit_p)
        events <- rbinom(1, arm_size, p)
        
        data_list[[length(data_list) + 1]] <- data.frame(
          study = paste0("Study_", i),
          treatment = study_treatments[j],
          n = arm_size,
          events = events,
          true_effect = true_effects[j]
        )
      } else {
        y <- study_baseline + true_effects[j] + rnorm(1, 0, 0.5)
        
        data_list[[length(data_list) + 1]] <- data.frame(
          study = paste0("Study_", i),
          treatment = study_treatments[j],
          n = arm_size,
          mean = y,
          sd = abs(rnorm(1, 1, 0.2)),
          true_effect = true_effects[j]
        )
      }
    }
  }
  
  data <- do.call(rbind, data_list)
  
  return(list(
    data = data,
    true_component_effects = component_effects,
    true_interaction_effects = interaction_effects,
    study_designs = study_designs
  ))
}

print("cNMA methodologies loaded successfully!")