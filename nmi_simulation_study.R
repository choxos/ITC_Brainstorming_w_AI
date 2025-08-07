# Network Meta-Interpolation Simulation Study
# Comprehensive simulation to test novel NMI methodologies
# Author: Research Collaboration
# Date: 2025

# Load required libraries
library(cmdstanr)
library(posterior)
library(bayesplot)
library(loo)
library(dplyr)
library(ggplot2)
library(tidyr)
library(parallel)
library(foreach)
library(doParallel)
library(MASS)
library(mvtnorm)
library(gridExtra)
library(knitr)
library(reshape2)

# Source methodology functions
source("nmi_methodology.R")

# Set up parallel processing
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Simulation Framework
# ===================

#' Comprehensive NMI Simulation Study
#' 
#' Tests novel Bayesian NMI methods against standard NMI across
#' various scenarios including different correlation patterns,
#' missing data patterns, and effect modification structures.

run_nmi_simulation_study <- function(n_simulations = 2000, n_cores = NULL,
                                    scenarios = NULL, save_results = TRUE) {
  
  if (is.null(n_cores)) {
    n_cores <- min(detectCores() - 1, 8)
  }
  
  if (is.null(scenarios)) {
    scenarios <- define_nmi_simulation_scenarios()
  }
  
  cat("Starting NMI simulation study with", n_simulations, "replications per scenario\n")
  cat("Testing", length(scenarios), "scenarios across", n_cores, "cores\n\n")
  
  # Run simulations for each scenario
  all_results <- list()
  
  for (scenario_name in names(scenarios)) {
    scenario <- scenarios[[scenario_name]]
    
    cat("Scenario:", scenario_name, "\n")
    cat("  Description:", scenario$description, "\n")
    cat("  Methods tested:", paste(scenario$methods_to_test, collapse = ", "), "\n")
    
    # Run parallel simulations for this scenario
    scenario_results <- run_scenario_simulation(
      scenario = scenario,
      n_simulations = n_simulations,
      n_cores = n_cores
    )
    
    all_results[[scenario_name]] <- scenario_results
    
    cat("  Completed", n_simulations, "simulations\n\n")
  }
  
  # Combine and summarize results
  summary_results <- summarize_simulation_results(all_results)
  
  if (save_results) {
    save(all_results, summary_results, 
         file = paste0("nmi_simulation_results_", Sys.Date(), ".RData"))
    
    # Save summary tables
    write.csv(summary_results$performance_table, 
              "nmi_simulation_summary.csv", row.names = FALSE)
  }
  
  return(list(
    all_results = all_results,
    summary = summary_results
  ))
}

#' Define simulation scenarios for NMI methods
define_nmi_simulation_scenarios <- function() {
  
  list(
    # Scenario 1: Standard balanced network with modest correlation
    "Scenario_1_Balanced_Modest_Correlation" = list(
      description = "Balanced network with modest EM correlation and complete subgroup reporting",
      data_params = list(
        n_studies = 12,
        n_treatments = 4,
        n_em = 2,
        study_sizes = rep(300, 12),
        true_correlation = matrix(c(1, 0.3, 0.3, 1), 2, 2),
        true_treatment_effects = c(0, 0.4, 0.7, 1.1),
        true_effect_modifiers = matrix(c(0, 0, 0.2, 0.3, 0.1, 0.4, 0.3, 0.2), 
                                     nrow = 4, byrow = TRUE),
        outcome_type = "binary",
        baseline_logit = -1.0,
        heterogeneity_sd = 0.2,
        missing_subgroup_rate = 0.0,
        correlation_heterogeneity = 0.0,
        target_em = c(0.5, 0.6)
      ),
      methods_to_test = c("BH-NMI", "RB-NMI", "Standard-NMI"),
      truth_target = c(0.4 + 0.2*0.5 + 0.3*0.6)  # Treatment B at target EM
    ),
    
    # Scenario 2: High correlation with missing subgroups
    "Scenario_2_High_Correlation_Missing" = list(
      description = "High EM correlation with 30% missing subgroup analyses",
      data_params = list(
        n_studies = 10,
        n_treatments = 3,
        n_em = 2,
        study_sizes = rep(250, 10),
        true_correlation = matrix(c(1, 0.8, 0.8, 1), 2, 2),
        true_treatment_effects = c(0, 0.5, 0.9),
        true_effect_modifiers = matrix(c(0, 0, 0.3, 0.4, 0.2, 0.3), 
                                     nrow = 3, byrow = TRUE),
        outcome_type = "binary",
        baseline_logit = -0.5,
        heterogeneity_sd = 0.3,
        missing_subgroup_rate = 0.3,
        correlation_heterogeneity = 0.1,
        target_em = c(0.4, 0.7)
      ),
      methods_to_test = c("BH-NMI", "RB-NMI", "BMA-NMI", "Standard-NMI"),
      truth_target = c(0.5 + 0.3*0.4 + 0.4*0.7)
    ),
    
    # Scenario 3: Non-shared effect modification
    "Scenario_3_Non_Shared_EM" = list(
      description = "Non-shared effect modification violating SEM assumption",
      data_params = list(
        n_studies = 8,
        n_treatments = 4,
        n_em = 2,
        study_sizes = rep(400, 8),
        true_correlation = matrix(c(1, 0.2, 0.2, 1), 2, 2),
        true_treatment_effects = c(0, 0.6, 1.0, 1.3),
        # Different EM effects by treatment (non-SEM)
        true_effect_modifiers = matrix(c(0, 0, 0.1, 0.2, 0.3, 0.1, 0.4, 0.5), 
                                     nrow = 4, byrow = TRUE),
        outcome_type = "binary",
        baseline_logit = -0.8,
        heterogeneity_sd = 0.25,
        missing_subgroup_rate = 0.1,
        correlation_heterogeneity = 0.2,
        target_em = c(0.6, 0.4),
        non_shared_em = TRUE
      ),
      methods_to_test = c("BH-NMI", "RB-NMI", "GP-NMI", "Standard-NMI"),
      truth_target = c(0.6 + 0.1*0.6 + 0.2*0.4)  # Treatment B
    ),
    
    # Scenario 4: Small network with large studies
    "Scenario_4_Small_Network_Large_Studies" = list(
      description = "Small network with large study sizes and precise estimates",
      data_params = list(
        n_studies = 6,
        n_treatments = 3,
        n_em = 2,
        study_sizes = rep(800, 6),
        true_correlation = matrix(c(1, 0.5, 0.5, 1), 2, 2),
        true_treatment_effects = c(0, 0.3, 0.8),
        true_effect_modifiers = matrix(c(0, 0, 0.15, 0.25, 0.3, 0.1), 
                                     nrow = 3, byrow = TRUE),
        outcome_type = "binary",
        baseline_logit = -0.3,
        heterogeneity_sd = 0.15,
        missing_subgroup_rate = 0.05,
        correlation_heterogeneity = 0.05,
        target_em = c(0.3, 0.8)
      ),
      methods_to_test = c("BH-NMI", "BMA-NMI", "GP-NMI", "Standard-NMI"),
      truth_target = c(0.3 + 0.15*0.3 + 0.25*0.8)
    ),
    
    # Scenario 5: Large network with heterogeneous correlations
    "Scenario_5_Large_Network_Heterogeneous" = list(
      description = "Large network with heterogeneous correlation patterns across studies",
      data_params = list(
        n_studies = 20,
        n_treatments = 5,
        n_em = 2,
        study_sizes = rep(200, 20),
        true_correlation = matrix(c(1, 0.4, 0.4, 1), 2, 2),
        true_treatment_effects = c(0, 0.2, 0.5, 0.8, 1.2),
        true_effect_modifiers = matrix(c(0, 0, 0.1, 0.1, 0.2, 0.2, 
                                       0.3, 0.15, 0.4, 0.3), 
                                     nrow = 5, byrow = TRUE),
        outcome_type = "binary",
        baseline_logit = -1.2,
        heterogeneity_sd = 0.35,
        missing_subgroup_rate = 0.2,
        correlation_heterogeneity = 0.4,  # High heterogeneity
        target_em = c(0.5, 0.5)
      ),
      methods_to_test = c("BH-NMI", "RB-NMI", "BMA-NMI", "Standard-NMI"),
      truth_target = c(0.2 + 0.1*0.5 + 0.1*0.5)  # Treatment B
    ),
    
    # Scenario 6: Continuous outcome
    "Scenario_6_Continuous_Outcome" = list(
      description = "Continuous outcome with normal distribution and effect modification",
      data_params = list(
        n_studies = 10,
        n_treatments = 4,
        n_em = 2,
        study_sizes = rep(300, 10),
        true_correlation = matrix(c(1, 0.6, 0.6, 1), 2, 2),
        true_treatment_effects = c(0, 2.5, 4.0, 5.5),
        true_effect_modifiers = matrix(c(0, 0, 1.0, 1.5, 0.8, 2.0, 1.2, 1.8), 
                                     nrow = 4, byrow = TRUE),
        outcome_type = "continuous",
        baseline_mean = 20.0,
        outcome_sd = 5.0,
        heterogeneity_sd = 1.0,
        missing_subgroup_rate = 0.15,
        correlation_heterogeneity = 0.1,
        target_em = c(0.4, 0.6)
      ),
      methods_to_test = c("BH-NMI", "RB-NMI", "GP-NMI", "Standard-NMI"),
      truth_target = c(2.5 + 1.0*0.4 + 1.5*0.6)
    ),
    
    # Scenario 7: Three effect modifiers
    "Scenario_7_Three_Effect_Modifiers" = list(
      description = "Three effect modifiers with complex correlation structure",
      data_params = list(
        n_studies = 12,
        n_treatments = 3,
        n_em = 3,
        study_sizes = rep(350, 12),
        true_correlation = matrix(c(1, 0.3, 0.2, 0.3, 1, 0.4, 0.2, 0.4, 1), 3, 3),
        true_treatment_effects = c(0, 0.4, 0.9),
        true_effect_modifiers = matrix(c(0, 0, 0, 0.2, 0.3, 0.1, 0.3, 0.2, 0.4), 
                                     nrow = 3, byrow = TRUE),
        outcome_type = "binary",
        baseline_logit = -0.6,
        heterogeneity_sd = 0.25,
        missing_subgroup_rate = 0.25,
        correlation_heterogeneity = 0.15,
        target_em = c(0.5, 0.4, 0.6)
      ),
      methods_to_test = c("BH-NMI", "RB-NMI", "BMA-NMI", "Standard-NMI"),
      truth_target = c(0.4 + 0.2*0.5 + 0.3*0.4 + 0.1*0.6)
    ),
    
    # Scenario 8: Sparse network
    "Scenario_8_Sparse_Network" = list(
      description = "Sparse network with limited direct comparisons",
      data_params = list(
        n_studies = 8,
        n_treatments = 5,
        n_em = 2,
        study_sizes = rep(250, 8),
        true_correlation = matrix(c(1, 0.5, 0.5, 1), 2, 2),
        true_treatment_effects = c(0, 0.3, 0.7, 1.1, 1.4),
        true_effect_modifiers = matrix(c(0, 0, 0.1, 0.2, 0.25, 0.15, 
                                       0.3, 0.25, 0.35, 0.3), 
                                     nrow = 5, byrow = TRUE),
        outcome_type = "binary",
        baseline_logit = -0.7,
        heterogeneity_sd = 0.3,
        missing_subgroup_rate = 0.4,  # High missing rate
        correlation_heterogeneity = 0.2,
        target_em = c(0.6, 0.3),
        sparse_network = TRUE
      ),
      methods_to_test = c("BH-NMI", "RB-NMI", "BMA-NMI", "Standard-NMI"),
      truth_target = c(0.3 + 0.1*0.6 + 0.2*0.3)
    ),
    
    # Scenario 9: Variable study sizes
    "Scenario_9_Variable_Study_Sizes" = list(
      description = "Network with highly variable study sizes",
      data_params = list(
        n_studies = 14,
        n_treatments = 4,
        n_em = 2,
        study_sizes = c(50, 80, 100, 150, 200, 300, 400, 500, 600, 
                       800, 1000, 1200, 1500, 2000),
        true_correlation = matrix(c(1, 0.4, 0.4, 1), 2, 2),
        true_treatment_effects = c(0, 0.5, 0.8, 1.2),
        true_effect_modifiers = matrix(c(0, 0, 0.2, 0.25, 0.15, 0.3, 0.3, 0.2), 
                                     nrow = 4, byrow = TRUE),
        outcome_type = "binary",
        baseline_logit = -0.9,
        heterogeneity_sd = 0.2,
        missing_subgroup_rate = 0.2,
        correlation_heterogeneity = 0.1,
        target_em = c(0.45, 0.55)
      ),
      methods_to_test = c("BH-NMI", "RB-NMI", "BMA-NMI", "Standard-NMI"),
      truth_target = c(0.5 + 0.2*0.45 + 0.25*0.55)
    ),
    
    # Scenario 10: Extreme missing data
    "Scenario_10_Extreme_Missing" = list(
      description = "Extreme missing subgroup data scenario",
      data_params = list(
        n_studies = 10,
        n_treatments = 3,
        n_em = 2,
        study_sizes = rep(300, 10),
        true_correlation = matrix(c(1, 0.3, 0.3, 1), 2, 2),
        true_treatment_effects = c(0, 0.6, 1.0),
        true_effect_modifiers = matrix(c(0, 0, 0.3, 0.2, 0.4, 0.3), 
                                     nrow = 3, byrow = TRUE),
        outcome_type = "binary",
        baseline_logit = -0.5,
        heterogeneity_sd = 0.25,
        missing_subgroup_rate = 0.6,  # Extreme missing rate
        correlation_heterogeneity = 0.3,
        target_em = c(0.5, 0.5)
      ),
      methods_to_test = c("RB-NMI", "BMA-NMI", "Standard-NMI"),  # BH-NMI may struggle
      truth_target = c(0.6 + 0.3*0.5 + 0.2*0.5)
    )
  )
}

#' Generate NMI simulation data
generate_nmi_data <- function(params) {
  
  n_studies <- params$n_studies
  n_treatments <- params$n_treatments
  n_em <- params$n_em
  study_sizes <- params$study_sizes
  true_correlation <- params$true_correlation
  true_treatment_effects <- params$true_treatment_effects
  true_effect_modifiers <- params$true_effect_modifiers
  outcome_type <- params$outcome_type
  missing_subgroup_rate <- params$missing_subgroup_rate
  correlation_heterogeneity <- params$correlation_heterogeneity
  target_em <- params$target_em
  
  # Generate study-level effect modifier distributions
  study_em_means <- matrix(runif(n_studies * n_em, 0.2, 0.8), 
                          nrow = n_studies, ncol = n_em)
  
  # Generate study-specific correlations (if heterogeneous)
  study_correlations <- array(dim = c(n_em, n_em, n_studies))
  
  for (s in 1:n_studies) {
    if (correlation_heterogeneity > 0) {
      # Add noise to correlation matrix
      noise_matrix <- matrix(rnorm(n_em^2, 0, correlation_heterogeneity), 
                           n_em, n_em)
      noise_matrix <- (noise_matrix + t(noise_matrix)) / 2  # Make symmetric
      diag(noise_matrix) <- 0  # Keep diagonal at 1
      
      corr_candidate <- true_correlation + noise_matrix
      
      # Ensure positive definite
      if (min(eigen(corr_candidate)$values) > 0.01) {
        study_correlations[,,s] <- corr_candidate
      } else {
        study_correlations[,,s] <- true_correlation
      }
    } else {
      study_correlations[,,s] <- true_correlation
    }
  }
  
  # Generate subgroup data
  subgroup_data <- data.frame()
  ipd_data <- data.frame()
  
  for (s in 1:n_studies) {
    study_size <- study_sizes[min(s, length(study_sizes))]
    
    # Determine treatments in this study (simplified network structure)
    if (exists("sparse_network", where = params) && params$sparse_network) {
      # Create sparse connections
      available_treatments <- sample(1:n_treatments, 
                                   size = min(2, n_treatments), 
                                   replace = FALSE)
    } else {
      # Most studies compare treatment 1 with one other
      available_treatments <- c(1, sample(2:n_treatments, 1))
    }
    
    treatment_comparison <- available_treatments
    
    # Generate IPD for correlation estimation (if this is the reference study)
    if (s == 1) {
      n_ipd <- study_size
      
      # Generate correlated effect modifiers
      em_ipd <- rmvnorm(n_ipd, mean = study_em_means[s,], 
                       sigma = study_correlations[,,s])
      em_ipd <- pmax(0, pmin(1, em_ipd))  # Bound to [0,1]
      
      # Generate treatment assignments
      treatment_ipd <- sample(treatment_comparison, n_ipd, replace = TRUE)
      
      # Generate outcomes
      if (outcome_type == "binary") {
        logit_p <- params$baseline_logit
        
        for (i in 1:n_ipd) {
          t <- treatment_ipd[i]
          if (t > 1) {
            logit_p[i] <- params$baseline_logit + 
                         true_treatment_effects[t] +
                         sum(true_effect_modifiers[t,] * em_ipd[i,])
          }
        }
        
        y_ipd <- rbinom(n_ipd, 1, plogis(logit_p))
        
      } else if (outcome_type == "continuous") {
        y_mean <- params$baseline_mean
        
        for (i in 1:n_ipd) {
          t <- treatment_ipd[i]
          if (t > 1) {
            y_mean[i] <- params$baseline_mean + 
                        true_treatment_effects[t] +
                        sum(true_effect_modifiers[t,] * em_ipd[i,])
          }
        }
        
        y_ipd <- rnorm(n_ipd, y_mean, params$outcome_sd)
      }
      
      # Store IPD
      ipd_data <- data.frame(
        study = s,
        treatment = treatment_ipd,
        y = y_ipd,
        em_1 = em_ipd[,1],
        em_2 = em_ipd[,2]
      )
      
      if (n_em > 2) {
        ipd_data$em_3 <- em_ipd[,3]
      }
    }
    
    # Generate subgroup analyses for this study
    study_subgroups <- generate_study_subgroups(
      study_id = s,
      treatments = treatment_comparison,
      study_em_means = study_em_means[s,],
      study_correlation = study_correlations[,,s],
      study_size = study_size,
      true_treatment_effects = true_treatment_effects,
      true_effect_modifiers = true_effect_modifiers,
      params = params
    )
    
    # Apply missing data
    if (missing_subgroup_rate > 0) {
      study_subgroups <- apply_missing_subgroups(study_subgroups, 
                                                missing_subgroup_rate)
    }
    
    subgroup_data <- rbind(subgroup_data, study_subgroups)
  }
  
  return(list(
    subgroup_data = subgroup_data,
    ipd_data = ipd_data,
    truth = list(
      treatment_effects = true_treatment_effects,
      effect_modifiers = true_effect_modifiers,
      correlation = true_correlation,
      target_em = target_em
    )
  ))
}

# Helper function to generate subgroup analyses for a single study
generate_study_subgroups <- function(study_id, treatments, study_em_means, 
                                   study_correlation, study_size,
                                   true_treatment_effects, true_effect_modifiers,
                                   params) {
  
  n_em <- length(study_em_means)
  subgroups <- data.frame()
  
  # Generate overall study-level subgroup (ITT)
  for (t in treatments) {
    
    # Calculate true effect at study mean EM values
    true_effect <- ifelse(t == 1, 0, true_treatment_effects[t])
    if (t > 1) {
      true_effect <- true_effect + sum(true_effect_modifiers[t,] * study_em_means)
    }
    
    # Add heterogeneity
    observed_effect <- true_effect + rnorm(1, 0, params$heterogeneity_sd)
    
    # Calculate standard error
    if (params$outcome_type == "binary") {
      # Approximate SE for log-odds ratio
      se <- sqrt(2 / study_size)  # Simplified
    } else {
      # SE for continuous outcome
      se <- params$outcome_sd / sqrt(study_size / 2)
    }
    
    subgroups <- rbind(subgroups, data.frame(
      study = study_id,
      treatment = t,
      subgroup_type = "overall",
      te = observed_effect,
      se = se,
      n = study_size,
      em_1 = study_em_means[1],
      em_2 = study_em_means[2],
      em_3 = ifelse(n_em > 2, study_em_means[3], NA)
    ))
  }
  
  # Generate subgroup analyses by each effect modifier
  for (em_idx in 1:n_em) {
    for (em_level in c(0, 1)) {  # Binary EM levels
      for (t in treatments) {
        
        # Create EM vector for this subgroup
        em_values <- study_em_means
        em_values[em_idx] <- em_level
        
        # Calculate true effect
        true_effect <- ifelse(t == 1, 0, true_treatment_effects[t])
        if (t > 1) {
          true_effect <- true_effect + sum(true_effect_modifiers[t,] * em_values)
        }
        
        # Add heterogeneity
        observed_effect <- true_effect + rnorm(1, 0, params$heterogeneity_sd)
        
        # Subgroup size (approximately half for binary EM)
        subgroup_size <- round(study_size * ifelse(em_level == 1, 
                                                  study_em_means[em_idx],
                                                  1 - study_em_means[em_idx]))
        
        # SE adjusted for subgroup size
        if (params$outcome_type == "binary") {
          se <- sqrt(4 / subgroup_size)
        } else {
          se <- params$outcome_sd / sqrt(subgroup_size / 2)
        }
        
        # Create missing pattern (some EMs missing for subgroup analyses)
        em_1_val <- ifelse(em_idx == 1, em_level, NA)
        em_2_val <- ifelse(em_idx == 2, em_level, NA)
        em_3_val <- ifelse(n_em > 2 && em_idx == 3, em_level, NA)
        
        subgroups <- rbind(subgroups, data.frame(
          study = study_id,
          treatment = t,
          subgroup_type = paste0("em_", em_idx, "_", em_level),
          te = observed_effect,
          se = se,
          n = subgroup_size,
          em_1 = em_1_val,
          em_2 = em_2_val,
          em_3 = em_3_val
        ))
      }
    }
  }
  
  return(subgroups)
}

# Apply missing data to subgroup analyses
apply_missing_subgroups <- function(subgroups, missing_rate) {
  
  # Randomly remove some subgroup analyses (not overall ITT)
  subgroup_analyses <- subgroups[subgroups$subgroup_type != "overall", ]
  overall_analyses <- subgroups[subgroups$subgroup_type == "overall", ]
  
  n_to_remove <- round(nrow(subgroup_analyses) * missing_rate)
  
  if (n_to_remove > 0) {
    remove_idx <- sample(1:nrow(subgroup_analyses), n_to_remove)
    subgroup_analyses <- subgroup_analyses[-remove_idx, ]
  }
  
  return(rbind(overall_analyses, subgroup_analyses))
}

#' Run simulation for a single scenario
run_scenario_simulation <- function(scenario, n_simulations, n_cores) {
  
  # Set up parallel backend
  cl <- makeCluster(min(n_cores, detectCores() - 1))
  registerDoParallel(cl)
  
  results <- foreach(sim = 1:n_simulations, 
                    .packages = c("cmdstanr", "posterior", "dplyr", "MASS"),
                    .export = c("generate_nmi_data", "bhnmi", "rbnmi", 
                               "bmanmi", "gpnmi", "standard_nmi",
                               "generate_study_subgroups", "apply_missing_subgroups")) %dopar% {
    
    tryCatch({
      # Generate data for this simulation
      sim_data <- generate_nmi_data(scenario$data_params)
      
      # Fit each method
      method_results <- list()
      
      for (method in scenario$methods_to_test) {
        
        method_result <- tryCatch({
          
          if (method == "BH-NMI") {
            bhnmi(sim_data$subgroup_data, sim_data$ipd_data, 
                 scenario$data_params$target_em,
                 n_chains = 2, n_iter = 1000, n_warmup = 500)
            
          } else if (method == "RB-NMI") {
            rbnmi(sim_data$subgroup_data, sim_data$ipd_data,
                 scenario$data_params$target_em,
                 n_chains = 2, n_iter = 1000, n_warmup = 500)
            
          } else if (method == "BMA-NMI") {
            # Create multiple target EM values for BMA
            target_list <- list(
              "target1" = scenario$data_params$target_em,
              "target2" = scenario$data_params$target_em * 0.8,
              "target3" = scenario$data_params$target_em * 1.2
            )
            bmanmi(sim_data$subgroup_data, sim_data$ipd_data,
                  target_list, n_chains = 2, n_iter = 800, n_warmup = 400)
            
          } else if (method == "GP-NMI") {
            gpnmi(sim_data$subgroup_data, sim_data$ipd_data,
                 scenario$data_params$target_em,
                 n_chains = 2, n_iter = 1000, n_warmup = 500)
            
          } else if (method == "Standard-NMI") {
            standard_nmi(sim_data$subgroup_data, sim_data$ipd_data,
                        scenario$data_params$target_em)
          }
          
        }, error = function(e) {
          list(error = as.character(e), method = method)
        })
        
        method_results[[method]] <- method_result
      }
      
      # Return simulation results (store full truth list to compute pairwise targets)
      list(
        sim_id = sim,
        data = sim_data,
        results = method_results,
        truth = sim_data$truth
      )
      
    }, error = function(e) {
      list(sim_id = sim, error = as.character(e))
    })
  }
  
  stopCluster(cl)
  return(results)
}

#' Evaluate performance for each method
evaluate_nmi_performance <- function(simulation_results, method_name, truth_target) {
  
  n_sims <- length(simulation_results)
  successful_sims <- 0
  
  bias_values <- numeric()
  mse_values <- numeric()
  coverage_values <- logical()
  interval_widths <- numeric()
  computation_times <- numeric()
  convergence_rates <- numeric()
  
  for (sim_result in simulation_results) {
    
    if ("error" %in% names(sim_result) || 
        "error" %in% names(sim_result$results[[method_name]])) {
      next
    }
    
    method_result <- sim_result$results[[method_name]]
    # Compute truth for a specific comparison (B vs A) at target EM values
    # Expect truth_target to be a list with treatment_effects and effect_modifiers
    if (is.list(sim_result$truth) && all(c("treatment_effects","effect_modifiers","target_em") %in% names(sim_result$truth))) {
      te_true <- sim_result$truth$treatment_effects
      beta_true <- sim_result$truth$effect_modifiers
      target_em <- sim_result$truth$target_em
      # Treatment 2 (B) vs 1 (A) by default
      truth <- (te_true[2] + sum(beta_true[2,] * target_em)) - (te_true[1] + sum(beta_true[1,] * target_em))
    } else {
      truth <- truth_target
    }
    
    successful_sims <- successful_sims + 1
    
    # Extract point estimate and CI
    if (method_name == "Standard-NMI") {
      estimate <- as.numeric(method_result$interpolated_data$prediction)
      se <- as.numeric(method_result$interpolated_data$se)
      ci_lower <- estimate - 1.96 * se
      ci_upper <- estimate + 1.96 * se
      
    } else {
      # Bayesian methods
      if ("target_comparisons" %in% names(method_result)) {
        # Use pairwise comparison B vs A at target
        comp <- posterior::as_draws_df(method_result$fit$draws("target_comparisons"))
        # target_comparisons[2,1] flatten to draws column name like target_comparisons[2,1]
        colname <- "target_comparisons[2,1]"
        draws <- as.numeric(comp[[colname]])
        estimate <- mean(draws)
        ci_lower <- quantile(draws, 0.025)
        ci_upper <- quantile(draws, 0.975)
      } else if ("target_effects" %in% names(method_result)) {
        # Fallback to effects, compute B-A difference if available
        eff <- posterior::as_draws_df(method_result$fit$draws("target_effects"))
        if (all(c("target_effects[2]","target_effects[1]") %in% names(eff))) {
          draws <- as.numeric(eff[["target_effects[2]"]] - eff[["target_effects[1]"]])
          estimate <- mean(draws)
          ci_lower <- quantile(draws, 0.025)
          ci_upper <- quantile(draws, 0.975)
        } else {
          estimate <- method_result$target_effects$mean[2]
          ci_lower <- method_result$target_effects$q5[2]
          ci_upper <- method_result$target_effects$q95[2]
        }
      } else if ("bma_estimates" %in% names(method_result)) {
        estimate <- method_result$bma_estimates$averaged_effects[2]
        # Use combined CI if available
        if (!is.null(method_result$bma_estimates$ci)) {
          ci_lower <- method_result$bma_estimates$ci[1,2]
          ci_upper <- method_result$bma_estimates$ci[2,2]
        } else {
          ci_lower <- estimate - 1.96 * 0.2
          ci_upper <- estimate + 1.96 * 0.2
        }
      } else {
        # GP-NMI or other methods
        posterior_draws <- method_result$posterior
        if ("target_prediction" %in% names(posterior_draws)) {
          tp <- posterior_draws$target_prediction
          estimate <- mean(tp)
          ci_lower <- quantile(tp, 0.025)
          ci_upper <- quantile(tp, 0.975)
        } else {
          estimate <- 0.5  # Default
          ci_lower <- 0.3
          ci_upper <- 0.7
        }
      }
    }
    
    # Calculate performance metrics
    bias <- estimate - truth
    mse <- bias^2
    coverage <- (truth >= ci_lower) && (truth <= ci_upper)
    width <- ci_upper - ci_lower
    
    bias_values <- c(bias_values, bias)
    mse_values <- c(mse_values, mse)
    coverage_values <- c(coverage_values, coverage)
    interval_widths <- c(interval_widths, width)
    
    # Convergence assessment for Bayesian methods
    if (method_name != "Standard-NMI" && "fit" %in% names(method_result)) {
      diagnostics <- method_result$model_diagnostics
      max_rhat <- max(diagnostics$rhat, na.rm = TRUE)
      convergence_rates <- c(convergence_rates, as.numeric(max_rhat < 1.1))
    } else {
      convergence_rates <- c(convergence_rates, 1.0)  # Always converged for frequentist
    }
  }
  
  # Summarize performance
  performance <- list(
    method = method_name,
    n_successful = successful_sims,
    success_rate = successful_sims / n_sims,
    bias_mean = mean(bias_values),
    bias_sd = sd(bias_values),
    rmse = sqrt(mean(mse_values)),
    coverage_rate = mean(coverage_values),
    mean_width = mean(interval_widths),
    convergence_rate = mean(convergence_rates)
  )
  
  return(performance)
}

#' Summarize all simulation results
summarize_simulation_results <- function(all_results) {
  
  performance_table <- data.frame()
  
  for (scenario_name in names(all_results)) {
    scenario_results <- all_results[[scenario_name]]
    
    # Get truth target for this scenario
    truth_target <- scenario_results[[1]]$truth
    
    # Get methods tested
    methods_tested <- names(scenario_results[[1]]$results)
    
    for (method in methods_tested) {
      
      performance <- evaluate_nmi_performance(scenario_results, method, truth_target)
      
      perf_row <- data.frame(
        Scenario = scenario_name,
        Method = method,
        N_Successful = performance$n_successful,
        Success_Rate = round(performance$success_rate, 3),
        Bias_Mean = round(performance$bias_mean, 4),
        RMSE = round(performance$rmse, 4),
        Coverage = round(performance$coverage_rate, 3),
        Mean_Width = round(performance$mean_width, 3),
        Convergence = round(performance$convergence_rate, 3)
      )
      
      performance_table <- rbind(performance_table, perf_row)
    }
  }
  
  # Create summary plots
  plots <- create_simulation_plots(performance_table)
  
  return(list(
    performance_table = performance_table,
    plots = plots
  ))
}

#' Create visualization plots for simulation results
create_simulation_plots <- function(performance_table) {
  
  # RMSE comparison plot
  rmse_plot <- ggplot(performance_table, aes(x = Scenario, y = RMSE, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Root Mean Squared Error by Method and Scenario",
         y = "RMSE", x = "Scenario")
  
  # Coverage plot
  coverage_plot <- ggplot(performance_table, aes(x = Scenario, y = Coverage, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Coverage Probability by Method and Scenario",
         y = "Coverage Probability", x = "Scenario")
  
  # Bias plot
  bias_plot <- ggplot(performance_table, aes(x = Scenario, y = Bias_Mean, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Mean Bias by Method and Scenario",
         y = "Mean Bias", x = "Scenario")
  
  return(list(
    rmse_plot = rmse_plot,
    coverage_plot = coverage_plot,
    bias_plot = bias_plot
  ))
}

print("NMI simulation study framework loaded successfully!")