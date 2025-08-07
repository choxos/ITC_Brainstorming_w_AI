# Demo Script for Advanced Bayesian MAIC Methods
# Author: Research Collaboration
# Date: 2025
# Purpose: Demonstrate novel MAIC methodologies with example data

# Load required libraries and methodology
source("maic_methodology.R")

library(ggplot2)
library(dplyr)
library(gridExtra)
library(knitr)

cat("=== Advanced Bayesian MAIC Methods Demo ===\n\n")

# Set seed for reproducibility
set.seed(12345)

#' Generate demonstration dataset
generate_demo_data <- function() {
  
  cat("Generating demonstration dataset...\n")
  
  # IPD Study (AB trial)
  n_ipd <- 300
  
  # Generate realistic covariates
  age <- rnorm(n_ipd, 65, 12)
  age <- pmax(18, pmin(90, age))  # Realistic age range
  
  sex <- rbinom(n_ipd, 1, 0.45)  # Slightly more males
  
  severity <- rnorm(n_ipd, 5.5, 2)
  severity <- pmax(0, pmin(10, severity))  # 0-10 severity scale
  
  comorbidity <- rbinom(n_ipd, 1, 0.3 + 0.01 * (age - 65))  # Age-related
  
  # Treatment assignment (1 = control, 2 = intervention)
  treatment <- sample(c(1, 2), n_ipd, replace = TRUE, prob = c(0.5, 0.5))
  
  # Generate outcomes with treatment effect and effect modification
  linear_pred <- -1.5 +  # Baseline logit
                0.6 * (treatment == 2) +  # Treatment effect
                0.02 * (age - 65) +  # Age effect
                0.3 * sex +  # Sex effect  
                0.1 * severity +  # Severity effect
                0.4 * comorbidity +  # Comorbidity effect
                0.015 * (treatment == 2) * (age - 65)  # Age-treatment interaction
  
  outcome <- rbinom(n_ipd, 1, plogis(linear_pred))
  
  ipd_data <- data.frame(
    treatment = treatment,
    outcome = outcome,
    age = age,
    sex = sex,
    severity = severity,
    comorbidity = comorbidity
  )
  
  # Target Population (different from IPD)
  target_population <- list(
    age_mean = 70,      # Older population
    sex_mean = 0.6,     # More males
    severity_mean = 6.2,  # More severe
    comorbidity_mean = 0.45,  # More comorbidities
    age_var = 150,
    sex_var = 0.24,
    severity_var = 4,
    comorbidity_var = 0.25
  )
  
  # AgD Study (AC trial) - for anchored comparison
  agd_data <- list(
    events = 89,  # Event count in control arm
    total_n = 250,  # Total patients
    age_mean = 68,
    sex_mean = 0.5,
    severity_mean = 5.8,
    comorbidity_mean = 0.35,
    age_var = 140,
    sex_var = 0.25,
    severity_var = 3.5,
    comorbidity_var = 0.23
  )
  
  # Multiple target populations for MT-MAIC demo
  target_populations <- list(
    # Young population
    list(
      age_mean = 55, sex_mean = 0.4, severity_mean = 4.5, comorbidity_mean = 0.2,
      age_var = 100, sex_var = 0.24, severity_var = 3, comorbidity_var = 0.16
    ),
    # Elderly population  
    list(
      age_mean = 75, sex_mean = 0.6, severity_mean = 7.0, comorbidity_mean = 0.6,
      age_var = 80, sex_var = 0.24, severity_var = 4, comorbidity_var = 0.24
    ),
    # Severe disease population
    list(
      age_mean = 67, sex_mean = 0.5, severity_mean = 8.5, comorbidity_mean = 0.5,
      age_var = 120, sex_var = 0.25, severity_var = 2, comorbidity_var = 0.25
    )
  )
  
  return(list(
    ipd_data = ipd_data,
    agd_data = agd_data,
    target_population = target_population,
    target_populations = target_populations
  ))
}

#' Demonstrate population overlap assessment
demonstrate_overlap <- function(demo_data) {
  
  cat("\n=== Population Overlap Assessment ===\n")
  
  ipd_cov <- demo_data$ipd_data[, c("age", "sex", "severity", "comorbidity")]
  target_means <- c(
    demo_data$target_population$age_mean,
    demo_data$target_population$sex_mean,
    demo_data$target_population$severity_mean,
    demo_data$target_population$comorbidity_mean
  )
  
  # Calculate Mahalanobis distances
  distances <- mahalanobis(ipd_cov, target_means, cov(ipd_cov))
  overlap_scores <- 1 / (1 + distances)
  
  cat("Population Overlap Summary:\n")
  cat("  Mean overlap score:", round(mean(overlap_scores), 3), "\n")
  cat("  Proportion with good overlap (>0.5):", 
      round(mean(overlap_scores > 0.5), 3), "\n")
  cat("  Extreme cases (<0.1):", sum(overlap_scores < 0.1), "patients\n")
  
  # Create overlap plot
  overlap_df <- data.frame(
    overlap_score = overlap_scores,
    treatment = factor(demo_data$ipd_data$treatment)
  )
  
  p_overlap <- ggplot(overlap_df, aes(x = overlap_score, fill = treatment)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
    labs(title = "Population Overlap Distribution",
         subtitle = "Red line indicates good overlap threshold",
         x = "Overlap Score", y = "Count") +
    theme_minimal()
  
  print(p_overlap)
  
  return(overlap_scores)
}

#' Run all MAIC methods on demo data
run_demo_analysis <- function(demo_data) {
  
  cat("\n=== Running MAIC Method Demonstrations ===\n")
  
  adjustment_vars <- c("age", "sex", "severity", "comorbidity")
  results <- list()
  
  # 1. Bayesian Hierarchical MAIC
  cat("\n1. Bayesian Hierarchical MAIC (BH-MAIC)...\n")
  
  tryCatch({
    bh_result <- bhmaic(
      ipd_data = demo_data$ipd_data,
      agd_data = demo_data$agd_data,
      target_population = demo_data$target_population,
      outcome_type = "binary",
      adjustment_vars = adjustment_vars,
      anchored = TRUE,
      hierarchical_strength = 1.0,
      n_chains = 2, n_iter = 1000, n_warmup = 500  # Reduced for demo
    )
    
    results$bh_maic <- bh_result
    cat("   ✓ BH-MAIC completed successfully\n")
    
    # Extract key results
    posterior_summary <- bh_result$target_comparisons
    target_comparison <- posterior_summary[posterior_summary$variable == "target_comparisons[2,1]", ]
    
    cat("   Treatment effect (log OR):", round(target_comparison$mean, 3), 
        "(", round(target_comparison$q5, 3), ",", round(target_comparison$q95, 3), ")\n")
    cat("   Effective sample size:", round(bh_result$effective_sample_size$mean, 0), "\n")
    
  }, error = function(e) {
    cat("   ✗ BH-MAIC failed:", e$message, "\n")
    results$bh_maic <- NULL
  })
  
  # 2. Multi-target MAIC
  cat("\n2. Multi-target MAIC (MT-MAIC)...\n")
  
  tryCatch({
    mt_result <- mtmaic(
      ipd_data = demo_data$ipd_data,
      agd_data = demo_data$agd_data,
      target_populations = demo_data$target_populations,
      outcome_type = "binary",
      adjustment_vars = adjustment_vars,
      anchored = TRUE,
      n_chains = 2, n_iter = 800, n_warmup = 400
    )
    
    results$mt_maic <- mt_result
    cat("   ✓ MT-MAIC completed successfully\n")
    cat("   Analyzed", mt_result$n_targets, "target populations simultaneously\n")
    
  }, error = function(e) {
    cat("   ✗ MT-MAIC failed:", e$message, "\n")
    results$mt_maic <- NULL
  })
  
  # 3. Robust MAIC
  cat("\n3. Robust MAIC (R-MAIC)...\n")
  
  tryCatch({
    r_result <- rmaic(
      ipd_data = demo_data$ipd_data,
      agd_data = demo_data$agd_data,
      target_population = demo_data$target_population,
      outcome_type = "binary",
      adjustment_vars = adjustment_vars,
      anchored = TRUE,
      robustness_method = "trimming",
      overlap_threshold = 0.1,
      n_chains = 2, n_iter = 800, n_warmup = 400
    )
    
    results$r_maic <- r_result
    cat("   ✓ R-MAIC completed successfully\n")
    cat("   Robustness method: trimming with threshold 0.1\n")
    
  }, error = function(e) {
    cat("   ✗ R-MAIC failed:", e$message, "\n")
    results$r_maic <- NULL
  })
  
  # 4. Standard MAIC for comparison
  cat("\n4. Standard MAIC (comparison)...\n")
  
  tryCatch({
    std_result <- standard_maic(
      ipd_data = demo_data$ipd_data,
      agd_data = demo_data$agd_data,
      target_population = demo_data$target_population,
      outcome_type = "binary",
      adjustment_vars = adjustment_vars,
      anchored = TRUE,
      robust_se = TRUE
    )
    
    results$standard_maic <- std_result
    cat("   ✓ Standard MAIC completed successfully\n")
    cat("   Treatment effect estimate:", round(std_result$comparisons$comparison, 3), "\n")
    cat("   Effective sample size:", round(std_result$effective_sample_size, 0), "\n")
    
  }, error = function(e) {
    cat("   ✗ Standard MAIC failed:", e$message, "\n")
    results$standard_maic <- NULL
  })
  
  return(results)
}

#' Compare methods and create visualizations
compare_methods <- function(results) {
  
  cat("\n=== Method Comparison ===\n")
  
  comparison_data <- list()
  
  # Extract estimates from each method
  if (!is.null(results$bh_maic)) {
    bh_summary <- results$bh_maic$target_comparisons
    bh_comparison <- bh_summary[bh_summary$variable == "target_comparisons[2,1]", ]
    
    comparison_data[["BH-MAIC"]] <- list(
      estimate = bh_comparison$mean,
      lower = bh_comparison$q5,
      upper = bh_comparison$q95,
      method = "Bayesian Hierarchical"
    )
  }
  
  if (!is.null(results$r_maic)) {
    r_summary <- results$r_maic$model_diagnostics
    r_comparison <- r_summary[grepl("target_comparisons\\[2,1\\]", r_summary$variable), ]
    
    if (nrow(r_comparison) > 0) {
      comparison_data[["R-MAIC"]] <- list(
        estimate = r_comparison$mean,
        lower = r_comparison$q5,
        upper = r_comparison$q95,
        method = "Robust"
      )
    }
  }
  
  if (!is.null(results$standard_maic)) {
    comparison_data[["Standard MAIC"]] <- list(
      estimate = results$standard_maic$comparisons$comparison,
      lower = results$standard_maic$comparisons$comparison - 1.96 * results$standard_maic$standard_errors[1],
      upper = results$standard_maic$comparisons$comparison + 1.96 * results$standard_maic$standard_errors[1],
      method = "Standard"
    )
  }
  
  if (length(comparison_data) > 0) {
    # Create comparison data frame
    comp_df <- do.call(rbind, lapply(names(comparison_data), function(name) {
      data.frame(
        method = name,
        estimate = comparison_data[[name]]$estimate,
        lower = comparison_data[[name]]$lower,
        upper = comparison_data[[name]]$upper,
        stringsAsFactors = FALSE
      )
    }))
    
    # Create forest plot
    p_forest <- ggplot(comp_df, aes(x = method, y = estimate)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
      coord_flip() +
      labs(title = "Treatment Effect Estimates by Method",
           subtitle = "Log odds ratio with 95% credible/confidence intervals",
           x = "Method", y = "Log Odds Ratio") +
      theme_minimal()
    
    print(p_forest)
    
    # Print comparison table
    cat("\nTreatment Effect Comparison:\n")
    comp_table <- comp_df
    comp_table[, 2:4] <- round(comp_table[, 2:4], 3)
    print(kable(comp_table, format = "simple"))
    
    return(comp_df)
  } else {
    cat("No valid results available for comparison\n")
    return(NULL)
  }
}

#' Generate diagnostic plots
create_diagnostics <- function(results) {
  
  cat("\n=== Diagnostic Plots ===\n")
  
  plots <- list()
  
  # BH-MAIC diagnostics
  if (!is.null(results$bh_maic) && "fit" %in% names(results$bh_maic)) {
    
    cat("Generating BH-MAIC diagnostics...\n")
    
    # Trace plots for key parameters
    tryCatch({
      p_trace <- bayesplot::mcmc_trace(results$bh_maic$fit$draws(c("alpha_mean", "beta_treatment")))
      plots$bh_trace <- p_trace
    }, error = function(e) {
      cat("Could not create trace plot:", e$message, "\n")
    })
    
    # R-hat plot
    tryCatch({
      rhat_values <- results$bh_maic$model_diagnostics$rhat
      rhat_df <- data.frame(
        parameter = seq_along(rhat_values),
        rhat = rhat_values
      )
      
      p_rhat <- ggplot(rhat_df, aes(x = parameter, y = rhat)) +
        geom_point() +
        geom_hline(yintercept = 1.01, color = "red", linetype = "dashed") +
        labs(title = "R-hat Convergence Diagnostics",
             subtitle = "Red line shows convergence threshold (1.01)",
             x = "Parameter Index", y = "R-hat") +
        theme_minimal()
      
      plots$rhat <- p_rhat
      
    }, error = function(e) {
      cat("Could not create R-hat plot:", e$message, "\n")
    })
  }
  
  # Weight distribution plots
  if (!is.null(results$standard_maic)) {
    
    cat("Generating weight distribution plot...\n")
    
    weights_df <- data.frame(
      weight = results$standard_maic$weights,
      treatment = factor(demo_data$ipd_data$treatment),
      patient_id = 1:length(results$standard_maic$weights)
    )
    
    p_weights <- ggplot(weights_df, aes(x = weight)) +
      geom_histogram(bins = 30, alpha = 0.7, fill = "steelblue") +
      geom_vline(xintercept = median(weights_df$weight), 
                color = "red", linetype = "dashed") +
      labs(title = "MAIC Weight Distribution",
           subtitle = "Red line shows median weight",
           x = "Weight", y = "Count") +
      theme_minimal()
    
    plots$weights <- p_weights
    print(p_weights)
  }
  
  return(plots)
}

#' Main demonstration function
main_demo <- function() {
  
  cat("Starting MAIC methods demonstration...\n")
  cat("This may take several minutes to complete.\n\n")
  
  # Generate demo data
  demo_data <<- generate_demo_data()  # Make global for other functions
  
  # Assess population overlap
  overlap_scores <- demonstrate_overlap(demo_data)
  
  # Run MAIC analyses
  results <- run_demo_analysis(demo_data)
  
  # Compare methods
  comparison <- compare_methods(results)
  
  # Create diagnostic plots
  diagnostics <- create_diagnostics(results)
  
  # Summary
  cat("\n=== Demo Summary ===\n")
  cat("Demonstration completed successfully!\n")
  cat("Key findings:\n")
  
  if (!is.null(results$bh_maic)) {
    cat("  • BH-MAIC provided full uncertainty quantification\n")
  }
  
  if (!is.null(results$mt_maic)) {
    cat("  • MT-MAIC analyzed", results$mt_maic$n_targets, "target populations\n")
  }
  
  if (!is.null(results$r_maic)) {
    cat("  • R-MAIC handled population overlap robustly\n")
  }
  
  cat("  • Population overlap score:", round(mean(overlap_scores), 3), "\n")
  
  if (!is.null(comparison)) {
    estimates_range <- range(comparison$estimate)
    cat("  • Treatment effect estimates ranged from", 
        round(estimates_range[1], 3), "to", round(estimates_range[2], 3), "\n")
  }
  
  cat("\nAll results and plots are available in the returned objects.\n")
  cat("For production use, increase the number of MCMC iterations.\n")
  
  return(list(
    demo_data = demo_data,
    results = results,
    comparison = comparison,
    diagnostics = diagnostics,
    overlap_scores = overlap_scores
  ))
}

# Run the demonstration
cat("Advanced Bayesian MAIC Methods - Interactive Demo\n")
cat("================================================\n\n")

demo_output <- main_demo()

cat("\n=== Demo Complete ===\n")
cat("Access results using: demo_output$results\n")
cat("View comparisons using: demo_output$comparison\n")
cat("See diagnostics using: demo_output$diagnostics\n")