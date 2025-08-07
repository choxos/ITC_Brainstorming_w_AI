# Comprehensive Simulation Study for Novel cNMA Methodologies
# Author: Research Collaboration
# Date: 2025
# Purpose: Test and compare new cNMA methods with extensive simulations

source("cnma_methodology.R")

# Simulation Study Framework
# =========================

#' Run comprehensive simulation study
#' 
#' @param n_simulations Number of simulation iterations
#' @param scenarios List of simulation scenarios
#' @param methods List of methods to compare
#' 
#' @return Results data frame

run_simulation_study <- function(n_simulations = 1000, scenarios = NULL, methods = NULL) {
  
  if (is.null(scenarios)) {
    scenarios <- define_simulation_scenarios()
  }
  
  if (is.null(methods)) {
    methods <- list(
      "BCI-NMA" = "bci_nma",
      "ACS-NMA" = "acs_nma", 
      "HCE-NMA" = "hce_nma",
      "Traditional" = "traditional_cnma"
    )
  }
  
  # Initialize results storage
  results <- list()
  
  # Run simulations for each scenario
  for (scenario_name in names(scenarios)) {
    cat("Running scenario:", scenario_name, "\n")
    scenario <- scenarios[[scenario_name]]
    
    scenario_results <- foreach(sim = 1:n_simulations, .combine = rbind,
                               .packages = c("cmdstanr", "posterior", "loo", "netmeta", "dplyr")) %dopar% {
      
      tryCatch({
        # Generate data for this simulation
        sim_data <- do.call(generate_cnma_data, scenario$data_params)
        
        # Test each method
        method_results <- list()
        
        for (method_name in names(methods)) {
          start_time <- Sys.time()
          
          if (method_name == "BCI-NMA") {
            result <- bci_nma(sim_data$data, scenario$component_matrix,
                             outcome_type = scenario$outcome_type,
                             n_chains = 2, n_iter = 1000, n_warmup = 500)
            
            performance <- evaluate_bci_nma_performance(result, sim_data)
            
          } else if (method_name == "ACS-NMA") {
            # Create alternative component matrices for model selection
            alt_matrices <- generate_alternative_matrices(scenario$component_matrix)
            
            result <- acs_nma(sim_data$data, alt_matrices,
                             outcome_type = scenario$outcome_type)
            
            performance <- evaluate_acs_nma_performance(result, sim_data, scenario)
            
          } else if (method_name == "HCE-NMA") {
            # Generate study covariates
            study_covariates <- matrix(rnorm(scenario$data_params$n_studies * 2), 
                                     ncol = 2)
            
            result <- hce_nma(sim_data$data, scenario$component_matrix,
                             study_covariates = study_covariates,
                             outcome_type = scenario$outcome_type)
            
            performance <- evaluate_hce_nma_performance(result, sim_data)
            
          } else if (method_name == "Traditional") {
            result <- traditional_cnma(sim_data$data, scenario$component_matrix,
                                     outcome_type = scenario$outcome_type)
            
            performance <- evaluate_traditional_performance(result, sim_data)
          }
          
          end_time <- Sys.time()
          
          method_results[[method_name]] <- data.frame(
            scenario = scenario_name,
            simulation = sim,
            method = method_name,
            bias = performance$bias,
            mse = performance$mse,
            coverage = performance$coverage,
            power = performance$power,
            computation_time = as.numeric(end_time - start_time),
            convergence = performance$convergence,
            model_selection_accuracy = performance$model_selection_accuracy %||% NA
          )
        }
        
        do.call(rbind, method_results)
        
      }, error = function(e) {
        # Return error row if simulation fails
        data.frame(
          scenario = scenario_name,
          simulation = sim,
          method = "ERROR",
          bias = NA, mse = NA, coverage = NA, power = NA,
          computation_time = NA, convergence = FALSE,
          model_selection_accuracy = NA
        )
      })
    }
    
    results[[scenario_name]] <- scenario_results
  }
  
  # Combine all results
  final_results <- do.call(rbind, results)
  
  return(final_results)
}

# Define Simulation Scenarios
# ===========================

define_simulation_scenarios <- function() {
  
  # Scenario 1: Simple additive effects (traditional cNMA should work well)
  component_matrix_1 <- matrix(c(
    1, 0, 0,  # Treatment A: Component 1 only
    0, 1, 0,  # Treatment B: Component 2 only
    0, 0, 1,  # Treatment C: Component 3 only
    1, 1, 0,  # Treatment D: Components 1+2
    1, 0, 1,  # Treatment E: Components 1+3
    0, 1, 1,  # Treatment F: Components 2+3
    1, 1, 1   # Treatment G: All components
  ), nrow = 7, ncol = 3, byrow = TRUE)
  
  rownames(component_matrix_1) <- paste0("Treatment_", LETTERS[1:7])
  colnames(component_matrix_1) <- paste0("Component_", 1:3)
  
  scenario_1 <- list(
    name = "Simple_Additive",
    component_matrix = component_matrix_1,
    data_params = list(
      n_studies = 20,
      treatments = rownames(component_matrix_1),
      component_matrix = component_matrix_1,
      component_effects = c(0.5, 0.3, 0.7),
      interaction_effects = NULL,
      outcome_type = "binary"
    ),
    outcome_type = "binary"
  )
  
  # Scenario 2: Strong synergistic interactions
  scenario_2 <- list(
    name = "Synergistic_Interactions",
    component_matrix = component_matrix_1,
    data_params = list(
      n_studies = 25,
      treatments = rownames(component_matrix_1),
      component_matrix = component_matrix_1,
      component_effects = c(0.2, 0.2, 0.2),
      interaction_effects = list(
        "1_2" = 0.8,  # Strong synergy between components 1 and 2
        "1_3" = 0.6,  # Moderate synergy between components 1 and 3
        "2_3" = 0.4   # Weak synergy between components 2 and 3
      ),
      outcome_type = "binary"
    ),
    outcome_type = "binary"
  )
  
  # Scenario 3: Antagonistic interactions
  scenario_3 <- list(
    name = "Antagonistic_Interactions", 
    component_matrix = component_matrix_1,
    data_params = list(
      n_studies = 25,
      treatments = rownames(component_matrix_1),
      component_matrix = component_matrix_1,
      component_effects = c(0.8, 0.8, 0.8),
      interaction_effects = list(
        "1_2" = -0.5,  # Antagonism between components 1 and 2
        "1_3" = -0.3,  # Weak antagonism between components 1 and 3
        "2_3" = -0.7   # Strong antagonism between components 2 and 3
      ),
      outcome_type = "binary"
    ),
    outcome_type = "binary"
  )
  
  # Scenario 4: Mixed interactions
  scenario_4 <- list(
    name = "Mixed_Interactions",
    component_matrix = component_matrix_1,
    data_params = list(
      n_studies = 30,
      treatments = rownames(component_matrix_1),
      component_matrix = component_matrix_1,
      component_effects = c(0.4, 0.4, 0.4),
      interaction_effects = list(
        "1_2" = 0.6,   # Synergy
        "1_3" = -0.3,  # Antagonism
        "2_3" = 0.0    # No interaction
      ),
      outcome_type = "binary"
    ),
    outcome_type = "binary"
  )
  
  # Scenario 5: Large network with many components
  component_matrix_5 <- matrix(0, nrow = 15, ncol = 6)
  # Create realistic combination patterns
  for (i in 1:6) {
    component_matrix_5[i, i] <- 1  # Single components
  }
  # Two-component combinations
  component_matrix_5[7, c(1,2)] <- 1
  component_matrix_5[8, c(1,3)] <- 1
  component_matrix_5[9, c(2,3)] <- 1
  component_matrix_5[10, c(4,5)] <- 1
  component_matrix_5[11, c(4,6)] <- 1
  component_matrix_5[12, c(5,6)] <- 1
  # Three-component combinations
  component_matrix_5[13, c(1,2,3)] <- 1
  component_matrix_5[14, c(4,5,6)] <- 1
  # All components
  component_matrix_5[15, ] <- 1
  
  rownames(component_matrix_5) <- paste0("Treatment_", 1:15)
  colnames(component_matrix_5) <- paste0("Component_", 1:6)
  
  scenario_5 <- list(
    name = "Large_Network",
    component_matrix = component_matrix_5,
    data_params = list(
      n_studies = 40,
      treatments = rownames(component_matrix_5),
      component_matrix = component_matrix_5,
      component_effects = c(0.3, 0.4, 0.2, 0.5, 0.3, 0.4),
      interaction_effects = list(
        "1_2" = 0.3,
        "2_3" = -0.2,
        "4_5" = 0.4,
        "5_6" = -0.1
      ),
      outcome_type = "binary"
    ),
    outcome_type = "binary"
  )
  
  # Scenario 6: Continuous outcome
  scenario_6 <- list(
    name = "Continuous_Outcome",
    component_matrix = component_matrix_1,
    data_params = list(
      n_studies = 20,
      treatments = rownames(component_matrix_1),
      component_matrix = component_matrix_1,
      component_effects = c(2.5, 1.8, 3.2),
      interaction_effects = list(
        "1_2" = 1.5,
        "1_3" = -0.8,
        "2_3" = 0.0
      ),
      outcome_type = "continuous"
    ),
    outcome_type = "continuous"
  )
  
  # Scenario 7: High heterogeneity
  scenario_7 <- list(
    name = "High_Heterogeneity",
    component_matrix = component_matrix_1,
    data_params = list(
      n_studies = 25,
      treatments = rownames(component_matrix_1),
      component_matrix = component_matrix_1,
      component_effects = c(0.5, 0.3, 0.7),
      interaction_effects = NULL,
      outcome_type = "binary"
    ),
    outcome_type = "binary",
    heterogeneity = "high"  # Will be used to add extra noise
  )
  
  # Scenario 8: Sparse network (disconnected components)
  component_matrix_8 <- matrix(c(
    1, 0, 0, 0,  # Group 1: Components 1
    1, 1, 0, 0,  # Group 1: Components 1+2
    0, 1, 0, 0,  # Group 1: Component 2
    0, 0, 1, 0,  # Group 2: Component 3  
    0, 0, 1, 1,  # Group 2: Components 3+4
    0, 0, 0, 1   # Group 2: Component 4
  ), nrow = 6, ncol = 4, byrow = TRUE)
  
  rownames(component_matrix_8) <- paste0("Treatment_", LETTERS[1:6])
  colnames(component_matrix_8) <- paste0("Component_", 1:4)
  
  scenario_8 <- list(
    name = "Sparse_Network",
    component_matrix = component_matrix_8,
    data_params = list(
      n_studies = 15,
      treatments = rownames(component_matrix_8),
      component_matrix = component_matrix_8,
      component_effects = c(0.4, 0.6, 0.3, 0.5),
      interaction_effects = NULL,
      outcome_type = "binary"
    ),
    outcome_type = "binary"
  )
  
  return(list(
    scenario_1 = scenario_1,
    scenario_2 = scenario_2,
    scenario_3 = scenario_3,
    scenario_4 = scenario_4,
    scenario_5 = scenario_5,
    scenario_6 = scenario_6,
    scenario_7 = scenario_7,
    scenario_8 = scenario_8
  ))
}

# Performance Evaluation Functions
# ================================

#' Evaluate BCI-NMA performance
evaluate_bci_nma_performance <- function(result, sim_data) {
  
  # Extract true and estimated component effects
  true_effects <- sim_data$true_component_effects
  estimated_effects <- result$component_effects$mean
  
  # Calculate bias
  bias <- mean(abs(estimated_effects - true_effects))
  
  # Calculate MSE
  mse <- mean((estimated_effects - true_effects)^2)
  
  # Calculate coverage
  coverage <- mean(
    (true_effects >= result$component_effects$q2.5) & 
    (true_effects <= result$component_effects$q97.5)
  )
  
  # Calculate power (proportion of truly non-zero effects detected)
  power <- mean(
    (true_effects != 0) & 
    ((result$component_effects$q2.5 > 0) | (result$component_effects$q97.5 < 0))
  )
  
  # Check convergence
  rhat_values <- summary(result$fit)$summary[, "Rhat"]
  convergence <- all(rhat_values < 1.1, na.rm = TRUE)
  
  return(list(
    bias = bias,
    mse = mse,
    coverage = coverage,
    power = power,
    convergence = convergence
  ))
}

#' Evaluate ACS-NMA performance
evaluate_acs_nma_performance <- function(result, sim_data, scenario) {
  
  # Model selection accuracy (did it choose the correct component structure?)
  best_model <- result$best_model
  
  # Evaluate the best model
  performance <- evaluate_bci_nma_performance(best_model, sim_data)
  
  # Add model selection accuracy
  # Simplified: assume correct model is the one with best WAIC
  model_selection_accuracy <- result$best_model_idx == 1  # Assume first model is correct
  
  performance$model_selection_accuracy <- model_selection_accuracy
  
  return(performance)
}

#' Evaluate HCE-NMA performance
evaluate_hce_nma_performance <- function(result, sim_data) {
  
  # Similar to BCI-NMA but focus on hierarchical effects
  posterior <- result$posterior
  
  # Extract main component effects
  true_effects <- sim_data$true_component_effects
  estimated_effects <- apply(posterior$beta_0, 2, mean)
  
  bias <- mean(abs(estimated_effects - true_effects))
  mse <- mean((estimated_effects - true_effects)^2)
  
  # Coverage for main effects
  coverage <- mean(
    (true_effects >= apply(posterior$beta_0, 2, quantile, 0.025)) & 
    (true_effects <= apply(posterior$beta_0, 2, quantile, 0.975))
  )
  
  # Power
  power <- mean(
    (true_effects != 0) & 
    ((apply(posterior$beta_0, 2, quantile, 0.025) > 0) | 
     (apply(posterior$beta_0, 2, quantile, 0.975) < 0))
  )
  
  # Check convergence
  rhat_values <- summary(result$fit)$summary[, "Rhat"]
  convergence <- all(rhat_values < 1.1, na.rm = TRUE)
  
  return(list(
    bias = bias,
    mse = mse,
    coverage = coverage,
    power = power,
    convergence = convergence
  ))
}

#' Evaluate traditional cNMA performance
evaluate_traditional_performance <- function(result, sim_data) {
  
  # Extract component effects from netcomb result
  true_effects <- sim_data$true_component_effects
  
  if (inherits(result, "netcomb")) {
    estimated_effects <- result$Comp.random
    comp_se <- result$seComp.random
    
    # Calculate 95% CI
    lower_ci <- estimated_effects - 1.96 * comp_se
    upper_ci <- estimated_effects + 1.96 * comp_se
    
    bias <- mean(abs(estimated_effects - true_effects))
    mse <- mean((estimated_effects - true_effects)^2)
    
    coverage <- mean((true_effects >= lower_ci) & (true_effects <= upper_ci))
    
    power <- mean(
      (true_effects != 0) & 
      ((lower_ci > 0) | (upper_ci < 0))
    )
    
    convergence <- TRUE  # Assume frequentist methods always converge
    
  } else {
    # If traditional method failed
    bias <- NA
    mse <- NA
    coverage <- NA
    power <- NA
    convergence <- FALSE
  }
  
  return(list(
    bias = bias,
    mse = mse,
    coverage = coverage,
    power = power,
    convergence = convergence
  ))
}

# Utility Functions for Simulation
# ================================

#' Generate alternative component matrices for model selection
generate_alternative_matrices <- function(true_matrix) {
  
  n_treatments <- nrow(true_matrix)
  n_components <- ncol(true_matrix)
  
  alternatives <- list()
  
  # Original matrix
  alternatives[[1]] <- true_matrix
  
  # Simpler matrix (fewer components)
  if (n_components > 2) {
    simple_matrix <- true_matrix[, 1:(n_components-1)]
    alternatives[[2]] <- simple_matrix
  }
  
  # More complex matrix (additional component)
  complex_matrix <- cbind(true_matrix, sample(0:1, n_treatments, replace = TRUE))
  alternatives[[3]] <- complex_matrix
  
  # Random matrix with same dimensions
  random_matrix <- matrix(sample(0:1, n_treatments * n_components, replace = TRUE),
                         nrow = n_treatments, ncol = n_components)
  alternatives[[4]] <- random_matrix
  
  return(alternatives)
}

#' Helper function for missing values
`%||%` <- function(x, y) if (is.null(x)) y else x

# Main Simulation Execution
# =========================

#' Run the complete simulation study
run_complete_simulation <- function(n_simulations = 500) {
  
  cat("Starting comprehensive cNMA simulation study...\n")
  cat("Number of simulations per scenario:", n_simulations, "\n")
  
  # Run simulations
  results <- run_simulation_study(n_simulations = n_simulations)
  
  # Save results
  saveRDS(results, "cnma_simulation_results.rds")
  
  # Calculate summary statistics
  summary_stats <- results %>%
    filter(method != "ERROR") %>%
    group_by(scenario, method) %>%
    summarise(
      mean_bias = mean(bias, na.rm = TRUE),
      mean_mse = mean(mse, na.rm = TRUE),
      mean_coverage = mean(coverage, na.rm = TRUE),
      mean_power = mean(power, na.rm = TRUE),
      convergence_rate = mean(convergence, na.rm = TRUE),
      mean_computation_time = mean(computation_time, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create performance plots
  create_performance_plots(results, summary_stats)
  
  # Create results tables
  create_results_tables(summary_stats)
  
  return(list(
    raw_results = results,
    summary_stats = summary_stats
  ))
}

# Visualization Functions
# =======================

create_performance_plots <- function(results, summary_stats) {
  
  # Plot 1: Bias comparison
  p1 <- ggplot(summary_stats, aes(x = method, y = mean_bias, fill = method)) +
    geom_bar(stat = "identity") +
    facet_wrap(~scenario, scales = "free") +
    theme_minimal() +
    labs(title = "Mean Absolute Bias by Method and Scenario",
         y = "Mean Absolute Bias", x = "Method") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("bias_comparison.png", p1, width = 12, height = 8)
  
  # Plot 2: MSE comparison
  p2 <- ggplot(summary_stats, aes(x = method, y = mean_mse, fill = method)) +
    geom_bar(stat = "identity") +
    facet_wrap(~scenario, scales = "free") +
    theme_minimal() +
    labs(title = "Mean Squared Error by Method and Scenario",
         y = "Mean Squared Error", x = "Method") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("mse_comparison.png", p2, width = 12, height = 8)
  
  # Plot 3: Coverage comparison
  p3 <- ggplot(summary_stats, aes(x = method, y = mean_coverage, fill = method)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    facet_wrap(~scenario) +
    theme_minimal() +
    labs(title = "95% Confidence Interval Coverage by Method and Scenario",
         y = "Coverage Probability", x = "Method") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("coverage_comparison.png", p3, width = 12, height = 8)
  
  # Plot 4: Computation time comparison
  p4 <- ggplot(summary_stats, aes(x = method, y = mean_computation_time, fill = method)) +
    geom_bar(stat = "identity") +
    facet_wrap(~scenario, scales = "free") +
    theme_minimal() +
    labs(title = "Mean Computation Time by Method and Scenario",
         y = "Computation Time (seconds)", x = "Method") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("computation_time_comparison.png", p4, width = 12, height = 8)
}

create_results_tables <- function(summary_stats) {
  
  # Main results table
  main_table <- summary_stats %>%
    arrange(scenario, method) %>%
    mutate(
      mean_bias = round(mean_bias, 3),
      mean_mse = round(mean_mse, 3),
      mean_coverage = round(mean_coverage, 3),
      mean_power = round(mean_power, 3),
      convergence_rate = round(convergence_rate, 3),
      mean_computation_time = round(mean_computation_time, 2)
    )
  
  write.csv(main_table, "cnma_simulation_summary.csv", row.names = FALSE)
  
  # Create LaTeX table
  latex_table <- xtable(main_table, 
                       caption = "Simulation Study Results for Component Network Meta-Analysis Methods",
                       label = "tab:simulation_results")
  
  print(latex_table, file = "simulation_results_table.tex", 
        include.rownames = FALSE, booktabs = TRUE)
}

cat("Simulation study framework loaded successfully!\n")
cat("Run run_complete_simulation() to execute the full study.\n")