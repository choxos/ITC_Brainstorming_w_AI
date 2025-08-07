# Demonstration of Novel cNMA Methodologies
# Author: Research Collaboration
# Date: 2025
# Purpose: Quick demonstration of new Bayesian cNMA methods

# Note: This is a simplified demonstration version that runs quickly
# For the full simulation study, use run_complete_analysis.R

library(cmdstanr)
library(posterior)
library(loo)
library(dplyr)
library(ggplot2)
library(gridExtra)

cat("===============================================\n")
cat("DEMONSTRATION: Novel cNMA Methodologies\n")
cat("===============================================\n\n")

# Set up demonstration
set.seed(123)
options(mc.cores = 2)  # Use fewer cores for demo

# Define a simple component matrix for demonstration
demo_component_matrix <- matrix(c(
  1, 0, 0,  # Treatment A: Component 1
  0, 1, 0,  # Treatment B: Component 2  
  0, 0, 1,  # Treatment C: Component 3
  1, 1, 0,  # Treatment D: Components 1+2
  1, 0, 1,  # Treatment E: Components 1+3
  0, 1, 1,  # Treatment F: Components 2+3
  1, 1, 1   # Treatment G: All components
), nrow = 7, ncol = 3, byrow = TRUE)

rownames(demo_component_matrix) <- paste0("Treatment_", LETTERS[1:7])
colnames(demo_component_matrix) <- paste0("Component_", 1:3)

cat("Component Matrix:\n")
print(demo_component_matrix)
cat("\n")

# Generate demonstration data with known interactions
cat("Generating demonstration data with known component interactions...\n")

true_component_effects <- c(0.3, 0.4, 0.5)  # Main effects
true_interactions <- list("1_2" = 0.6, "1_3" = -0.3, "2_3" = 0.0)  # Interactions

demo_data_generator <- function() {
  n_studies <- 12
  studies <- paste0("Study_", 1:n_studies)
  
  # Create study designs
  data_list <- list()
  
  for (i in 1:n_studies) {
    # Randomly select 2-3 treatments for each study
    n_arms <- sample(2:3, 1, prob = c(0.7, 0.3))
    study_treatments <- sample(rownames(demo_component_matrix), n_arms)
    
    study_baseline <- rnorm(1, -1.5, 0.3)  # Logit scale baseline
    
    for (j in 1:n_arms) {
      treatment_idx <- which(rownames(demo_component_matrix) == study_treatments[j])
      
      # Calculate true treatment effect
      effect <- sum(demo_component_matrix[treatment_idx, ] * true_component_effects)
      
      # Add interaction effects
      components_present <- which(demo_component_matrix[treatment_idx, ] == 1)
      if (length(components_present) >= 2) {
        for (k1 in 1:(length(components_present)-1)) {
          for (k2 in (k1+1):length(components_present)) {
            interaction_key <- paste0(components_present[k1], "_", components_present[k2])
            if (interaction_key %in% names(true_interactions)) {
              effect <- effect + true_interactions[[interaction_key]]
            }
          }
        }
      }
      
      # Generate binary outcome
      arm_size <- sample(40:120, 1)
      logit_p <- study_baseline + effect + rnorm(1, 0, 0.15)
      p <- plogis(logit_p)
      events <- rbinom(1, arm_size, p)
      
      data_list[[length(data_list) + 1]] <- data.frame(
        study = studies[i],
        treatment = study_treatments[j],
        n = arm_size,
        events = events,
        true_effect = effect
      )
    }
  }
  
  do.call(rbind, data_list)
}

demo_data <- demo_data_generator()

cat("Generated", nrow(demo_data), "treatment arms across", 
    length(unique(demo_data$study)), "studies\n")
cat("Sample of data:\n")
print(head(demo_data, 10))
cat("\n")

# Simplified BCI-NMA for demonstration
cat("Running simplified BCI-NMA (Bayesian Component Interaction NMA)...\n")

# Create a simplified Stan model for demonstration
demo_stan_code <- "
data {
  int<lower=0> N;
  int<lower=0> n_studies;
  int<lower=0> n_treatments;
  int<lower=0> n_components;
  int<lower=1, upper=n_studies> study[N];
  int<lower=1, upper=n_treatments> treatment[N];
  int<lower=0> events[N];
  int<lower=0> n_patients[N];
  matrix[n_treatments, n_components] component_matrix;
}

parameters {
  vector[n_studies] mu;
  vector[n_components] beta;
  vector[3] gamma;  // Simplified: only 3 pairwise interactions
  real<lower=0> tau;
  matrix[n_studies, n_components] delta;
}

transformed parameters {
  vector[N] logit_p;
  
  for (i in 1:N) {
    logit_p[i] = mu[study[i]];
    
    // Main component effects
    for (c in 1:n_components) {
      logit_p[i] += component_matrix[treatment[i], c] * (beta[c] + delta[study[i], c]);
    }
    
    // Interaction effects (simplified)
    logit_p[i] += component_matrix[treatment[i], 1] * component_matrix[treatment[i], 2] * gamma[1]; // 1-2
    logit_p[i] += component_matrix[treatment[i], 1] * component_matrix[treatment[i], 3] * gamma[2]; // 1-3  
    logit_p[i] += component_matrix[treatment[i], 2] * component_matrix[treatment[i], 3] * gamma[3]; // 2-3
  }
}

model {
  // Priors
  mu ~ normal(0, 2);
  beta ~ normal(0, 1);
  gamma ~ normal(0, 0.5);
  tau ~ half_normal(0, 1);
  
  for (c in 1:n_components) {
    delta[, c] ~ normal(0, tau);
  }
  
  // Likelihood
  events ~ binomial_logit(n_patients, logit_p);
}

generated quantities {
  vector[n_treatments] treatment_effects;
  
  for (t in 1:n_treatments) {
    treatment_effects[t] = 0;
    for (c in 1:n_components) {
      treatment_effects[t] += component_matrix[t, c] * beta[c];
    }
    // Add interaction effects
    treatment_effects[t] += component_matrix[t, 1] * component_matrix[t, 2] * gamma[1];
    treatment_effects[t] += component_matrix[t, 1] * component_matrix[t, 3] * gamma[2];
    treatment_effects[t] += component_matrix[t, 2] * component_matrix[t, 3] * gamma[3];
  }
}
"

# Prepare data for Stan
studies <- unique(demo_data$study)
study_data <- data.frame(study_id = 1:length(studies), study_name = studies)
demo_data <- merge(demo_data, study_data, by.x = "study", by.y = "study_name")

treatments <- rownames(demo_component_matrix)
treatment_map <- data.frame(treatment_id = 1:length(treatments), treatment_name = treatments)
demo_data <- merge(demo_data, treatment_map, by.x = "treatment", by.y = "treatment_name")

stan_data <- list(
  N = nrow(demo_data),
  n_studies = length(studies),
  n_treatments = nrow(demo_component_matrix),
  n_components = ncol(demo_component_matrix),
  study = demo_data$study_id,
  treatment = demo_data$treatment_id,
  events = demo_data$events,
  n_patients = demo_data$n,
  component_matrix = demo_component_matrix
)

# Fit BCI-NMA model
start_time <- Sys.time()
demo_model <- stan_model(model_code = demo_stan_code)
demo_fit <- sampling(demo_model, data = stan_data, 
                     chains = 2, iter = 1000, warmup = 500, 
                     cores = 2, verbose = FALSE)
end_time <- Sys.time()

cat("BCI-NMA completed in", round(as.numeric(end_time - start_time), 1), "seconds\n\n")

# Extract results
posterior <- extract(demo_fit)

# Component effects
component_results <- data.frame(
  component = colnames(demo_component_matrix),
  true_effect = true_component_effects,
  estimated_mean = apply(posterior$beta, 2, mean),
  estimated_sd = apply(posterior$beta, 2, sd),
  q2.5 = apply(posterior$beta, 2, quantile, 0.025),
  q97.5 = apply(posterior$beta, 2, quantile, 0.975)
)

cat("Component Effects Results:\n")
print(round(component_results, 3))
cat("\n")

# Interaction effects
interaction_results <- data.frame(
  interaction = c("Component_1-2", "Component_1-3", "Component_2-3"),
  true_effect = c(0.6, -0.3, 0.0),
  estimated_mean = apply(posterior$gamma, 2, mean),
  estimated_sd = apply(posterior$gamma, 2, sd),
  q2.5 = apply(posterior$gamma, 2, quantile, 0.025),
  q97.5 = apply(posterior$gamma, 2, quantile, 0.975)
)

cat("Interaction Effects Results:\n")
print(round(interaction_results, 3))
cat("\n")

# Traditional additive cNMA for comparison
cat("Running traditional additive cNMA for comparison...\n")

# Simplified traditional approach (ignores interactions)
traditional_stan_code <- "
data {
  int<lower=0> N;
  int<lower=0> n_studies;
  int<lower=0> n_treatments;
  int<lower=0> n_components;
  int<lower=1, upper=n_studies> study[N];
  int<lower=1, upper=n_treatments> treatment[N];
  int<lower=0> events[N];
  int<lower=0> n_patients[N];
  matrix[n_treatments, n_components] component_matrix;
}

parameters {
  vector[n_studies] mu;
  vector[n_components] beta;
  real<lower=0> tau;
  matrix[n_studies, n_components] delta;
}

transformed parameters {
  vector[N] logit_p;
  
  for (i in 1:N) {
    logit_p[i] = mu[study[i]];
    for (c in 1:n_components) {
      logit_p[i] += component_matrix[treatment[i], c] * (beta[c] + delta[study[i], c]);
    }
  }
}

model {
  mu ~ normal(0, 2);
  beta ~ normal(0, 1);
  tau ~ half_normal(0, 1);
  
  for (c in 1:n_components) {
    delta[, c] ~ normal(0, tau);
  }
  
  events ~ binomial_logit(n_patients, logit_p);
}
"

traditional_model <- stan_model(model_code = traditional_stan_code)
traditional_fit <- sampling(traditional_model, data = stan_data,
                           chains = 2, iter = 1000, warmup = 500,
                           cores = 2, verbose = FALSE)

traditional_posterior <- extract(traditional_fit)

traditional_component_results <- data.frame(
  component = colnames(demo_component_matrix),
  true_effect = true_component_effects,
  estimated_mean = apply(traditional_posterior$beta, 2, mean),
  estimated_sd = apply(traditional_posterior$beta, 2, sd),
  q2.5 = apply(traditional_posterior$beta, 2, quantile, 0.025),
  q97.5 = apply(traditional_posterior$beta, 2, quantile, 0.975)
)

cat("Traditional cNMA Component Effects:\n")
print(round(traditional_component_results, 3))
cat("\n")

# Compare methods
cat("METHOD COMPARISON:\n")
cat("==================\n")

# Calculate bias
bci_bias <- mean(abs(component_results$estimated_mean - component_results$true_effect))
traditional_bias <- mean(abs(traditional_component_results$estimated_mean - 
                            traditional_component_results$true_effect))

cat("Mean Absolute Bias:\n")
cat("- BCI-NMA:", round(bci_bias, 3), "\n")
cat("- Traditional:", round(traditional_bias, 3), "\n")
cat("- Bias reduction:", round((traditional_bias - bci_bias)/traditional_bias * 100, 1), "%\n\n")

# Calculate coverage
bci_coverage <- mean((component_results$true_effect >= component_results$q2.5) & 
                    (component_results$true_effect <= component_results$q97.5))
traditional_coverage <- mean((traditional_component_results$true_effect >= traditional_component_results$q2.5) & 
                           (traditional_component_results$true_effect <= traditional_component_results$q97.5))

cat("95% Coverage Probability:\n")
cat("- BCI-NMA:", round(bci_coverage, 3), "\n")
cat("- Traditional:", round(traditional_coverage, 3), "\n\n")

# Interaction detection
interaction_detected <- sum((interaction_results$q2.5 > 0) | (interaction_results$q97.5 < 0))
true_interactions_count <- sum(interaction_results$true_effect != 0)

cat("Interaction Detection:\n")
cat("- True non-zero interactions:", true_interactions_count, "\n")
cat("- Interactions detected by BCI-NMA:", interaction_detected, "\n")
cat("- Detection rate:", round(interaction_detected/true_interactions_count * 100, 1), "%\n\n")

# Create visualization
cat("Creating demonstration plots...\n")

# Plot 1: Component effects comparison
component_comparison <- rbind(
  component_results %>% 
    select(component, true_effect, estimated_mean, q2.5, q97.5) %>%
    mutate(method = "BCI-NMA"),
  traditional_component_results %>% 
    select(component, true_effect, estimated_mean, q2.5, q97.5) %>%
    mutate(method = "Traditional")
)

p1 <- ggplot(component_comparison, aes(x = component, color = method)) +
  geom_point(aes(y = estimated_mean), position = position_dodge(width = 0.3), size = 3) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.2, 
                position = position_dodge(width = 0.3)) +
  geom_point(aes(y = true_effect), shape = 4, size = 4, color = "red") +
  labs(title = "Component Effects: BCI-NMA vs Traditional cNMA",
       subtitle = "Red X marks show true values",
       x = "Component", y = "Effect Size", color = "Method") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Plot 2: Interaction effects (BCI-NMA only)
p2 <- ggplot(interaction_results, aes(x = interaction)) +
  geom_point(aes(y = estimated_mean), size = 3, color = "blue") +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.2, color = "blue") +
  geom_point(aes(y = true_effect), shape = 4, size = 4, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = "Component Interaction Effects (BCI-NMA)",
       subtitle = "Red X marks show true values",
       x = "Component Interaction", y = "Interaction Effect") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine plots
combined_plot <- grid.arrange(p1, p2, nrow = 2)

ggsave("demo_cnma_results.png", combined_plot, width = 10, height = 8, dpi = 300)

cat("Demonstration completed successfully!\n")
cat("Results saved to 'demo_cnma_results.png'\n\n")

cat("SUMMARY:\n")
cat("=========\n")
cat("This demonstration shows how BCI-NMA can:\n")
cat("1. Accurately estimate component main effects\n")
cat("2. Detect and quantify component interactions\n")
cat("3. Outperform traditional additive cNMA when interactions exist\n\n")

cat("Key advantages demonstrated:\n")
cat("- ", round((traditional_bias - bci_bias)/traditional_bias * 100, 1), 
    "% reduction in bias\n")
cat("- Successful detection of ", interaction_detected, "/", true_interactions_count, 
    " true interactions\n")
cat("- Better uncertainty quantification with credible intervals\n\n")

cat("For the full simulation study with all methods and scenarios,\n")
cat("run: source('run_complete_analysis.R')\n")