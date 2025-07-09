# Regularized BMR Parameter Estimation with Coordinate Descent
# Hierarchical modeling with penalty on deviation from global parameters
# Y = a_gender + b_gender * W^0.75, with shrinkage toward joint estimates

# Load required libraries
library(data.table)
library(ggplot2)

# ============================================================================
# Main function returning all handles
# ============================================================================
bmr_regularized_analysis <- function() {
  handles <- list()
  handles$generate_hb_data <- generate_hb_data
  handles$fit_hierarchical_model <- fit_hierarchical_model
  handles$coordinate_descent <- coordinate_descent
  handles$regularization_path <- regularization_path
  handles$cross_validate <- cross_validate
  handles$plot_results <- plot_results
  handles$run_full_analysis <- run_full_analysis
  
  cat("Usage: h <- bmr_regularized_analysis()\n")
  cat("Quick start: results <- h$run_full_analysis()\n")
  cat("Interactive: data <- h$generate_hb_data(); fit <- h$fit_hierarchical_model(data)\n")
  
  return(handles)
}

# ============================================================================
# Data Generation from Harris-Benedict Equations
# ============================================================================
generate_hb_data <- function(weight_range = c(40, 100),    # kg
                            height_range = c(150, 200),    # cm  
                            age = 50,                      # years
                            n_points = 50) {
  
  # Create grid of weight and height values
  weights <- seq(weight_range[1], weight_range[2], length.out = n_points)
  heights <- seq(height_range[1], height_range[2], length.out = n_points)
  
  # Generate all combinations
  grid <- expand.grid(weight = weights, height = heights)
  
  # Harris-Benedict equations (original 1919 version)
  # Males: BMR = 66.473 + (13.7516 × weight_kg) + (5.0033 × height_cm) - (6.755 × age)
  # Females: BMR = 655.0955 + (9.5634 × weight_kg) + (1.8500 × height_cm) - (4.676 × age)
  
  # Calculate BMR for both genders
  bmr_male <- 66.473 + (13.7516 * grid$weight) + (5.0033 * grid$height) - (6.755 * age)
  bmr_female <- 655.0955 + (9.5634 * grid$weight) + (1.8500 * grid$height) - (4.676 * age)
  
  # Create data.table with weight^0.75 predictor
  dt_male <- data.table(
    gender = "Male",
    weight = grid$weight,           # kg
    height = grid$height,          # cm
    weight_075 = grid$weight^0.75, # kg^0.75 (dimensionless in this context)
    bmr = bmr_male                 # kcal/day
  )
  
  dt_female <- data.table(
    gender = "Female", 
    weight = grid$weight,
    height = grid$height,
    weight_075 = grid$weight^0.75,
    bmr = bmr_female
  )
  
  # Combine datasets
  dt_combined <- rbind(dt_male, dt_female)
  
  # Add some realistic measurement noise (CV ~5% for BMR measurements)
  set.seed(42)
  dt_combined[, bmr_noisy := bmr + rnorm(.N, mean = 0, sd = 0.05 * bmr)]
  
  cat("Generated", nrow(dt_combined), "data points\n")
  cat("Weight range:", range(dt_combined$weight), "kg\n")
  cat("BMR range:", round(range(dt_combined$bmr)), "kcal/day\n")
  
  return(dt_combined)
}

# ============================================================================
# Hierarchical Model Fitting with Coordinate Descent
# ============================================================================
fit_hierarchical_model <- function(dt, 
                                  lambda_a = 1.0,      # penalty weight for intercept
                                  lambda_b = 1.0,      # penalty weight for slope  
                                  max_iter = 100,
                                  tolerance = 1e-6,
                                  use_noisy = FALSE) {
  
  # Choose response variable
  y_var <- if(use_noisy) "bmr_noisy" else "bmr"
  
  # Fit global model first (both genders combined)
  global_fit <- lm(get(y_var) ~ weight_075, data = dt)
  a_global <- coef(global_fit)[1]  # intercept (kcal/day)
  b_global <- coef(global_fit)[2]  # slope (kcal/day per kg^0.75)
  
  cat("Global parameters:\n")
  cat("  a_global =", round(a_global, 2), "kcal/day\n")
  cat("  b_global =", round(b_global, 2), "kcal/day per kg^0.75\n")
  
  # Split data by gender
  dt_m <- dt[gender == "Male"]
  dt_f <- dt[gender == "Female"]
  
  # Design matrices
  X_m <- cbind(1, dt_m$weight_075)  # [1, W^0.75]
  X_f <- cbind(1, dt_f$weight_075)
  y_m <- dt_m[[y_var]]
  y_f <- dt_f[[y_var]]
  
  # Run coordinate descent
  result <- coordinate_descent(X_m, y_m, X_f, y_f, 
                              a_global, b_global,
                              lambda_a, lambda_b, 
                              max_iter, tolerance)
  
  # Calculate final predictions and residuals
  pred_m <- X_m %*% c(result$a_male, result$b_male)
  pred_f <- X_f %*% c(result$a_female, result$b_female)
  
  # Performance metrics
  rmse_m <- sqrt(mean((y_m - pred_m)^2))
  rmse_f <- sqrt(mean((y_f - pred_f)^2))
  rmse_overall <- sqrt(mean(c((y_m - pred_m)^2, (y_f - pred_f)^2)))
  
  r2_m <- 1 - sum((y_m - pred_m)^2) / sum((y_m - mean(y_m))^2)
  r2_f <- 1 - sum((y_f - pred_f)^2) / sum((y_f - mean(y_f))^2)
  
  # Create comprehensive results
  results <- list(
    # Parameters
    a_male = result$a_male,
    b_male = result$b_male,
    a_female = result$a_female,
    b_female = result$b_female,
    a_global = a_global,
    b_global = b_global,
    
    # Regularization
    lambda_a = lambda_a,
    lambda_b = lambda_b,
    
    # Convergence info
    converged = result$converged,
    iterations = result$iterations,
    convergence_path = result$path,
    
    # Performance
    rmse_male = rmse_m,
    rmse_female = rmse_f,
    rmse_overall = rmse_overall,
    r2_male = r2_m,
    r2_female = r2_f,
    
    # Data
    data = dt,
    predictions = data.table(
      gender = c(rep("Male", length(pred_m)), rep("Female", length(pred_f))),
      observed = c(y_m, y_f),
      predicted = c(pred_m, pred_f),
      weight_075 = c(dt_m$weight_075, dt_f$weight_075)
    )
  )
  
  class(results) <- "bmr_hierarchical_fit"
  return(results)
}

# ============================================================================
# Coordinate Descent Algorithm
# ============================================================================
coordinate_descent <- function(X_m, y_m, X_f, y_f,
                              a_global, b_global,
                              lambda_a, lambda_b,
                              max_iter = 100, 
                              tolerance = 1e-6) {
  
  # Initialize parameters at global values
  a_m <- a_global
  b_m <- b_global  
  a_f <- a_global
  b_f <- b_global
  
  # Storage for convergence monitoring
  path <- data.table(
    iteration = integer(),
    a_male = numeric(),
    b_male = numeric(), 
    a_female = numeric(),
    b_female = numeric(),
    objective = numeric()
  )
  
  for(iter in 1:max_iter) {
    # Store previous values
    a_m_old <- a_m; b_m_old <- b_m
    a_f_old <- a_f; b_f_old <- b_f
    
    # ========================================================================
    # Step 1: Update intercepts (a_m, a_f) with slopes fixed
    # ========================================================================
    
    # Residuals with current slopes
    resid_m_a <- y_m - b_m * X_m[,2]  # y - b*W^0.75
    resid_f_a <- y_f - b_f * X_f[,2]
    
    # Normal equations with regularization
    # For males: (n + lambda_a) * a_m = sum(resid_m_a) + lambda_a * a_global
    n_m <- length(y_m)
    n_f <- length(y_f)
    
    a_m <- (sum(resid_m_a) + lambda_a * a_global) / (n_m + lambda_a)
    a_f <- (sum(resid_f_a) + lambda_a * a_global) / (n_f + lambda_a)
    
    # ========================================================================  
    # Step 2: Update slopes (b_m, b_f) with intercepts fixed
    # ========================================================================
    
    # Residuals with current intercepts
    resid_m_b <- y_m - a_m  # y - a
    resid_f_b <- y_f - a_f
    
    # Normal equations for slopes
    # b_m: (sum(W^1.5) + lambda_b) * b_m = sum(W^0.75 * resid_m_b) + lambda_b * b_global
    sum_w075_sq_m <- sum(X_m[,2]^2)  # sum(W^1.5)
    sum_w075_sq_f <- sum(X_f[,2]^2)
    
    b_m <- (sum(X_m[,2] * resid_m_b) + lambda_b * b_global) / (sum_w075_sq_m + lambda_b)
    b_f <- (sum(X_f[,2] * resid_f_b) + lambda_b * b_global) / (sum_w075_sq_f + lambda_b)
    
    # ========================================================================
    # Calculate objective function
    # ========================================================================
    pred_m <- a_m + b_m * X_m[,2]
    pred_f <- a_f + b_f * X_f[,2]
    
    mse <- mean(c((y_m - pred_m)^2, (y_f - pred_f)^2))
    penalty <- lambda_a * ((a_m - a_global)^2 + (a_f - a_global)^2) + 
               lambda_b * ((b_m - b_global)^2 + (b_f - b_global)^2)
    objective <- mse + penalty
    
    # Store path
    path <- rbind(path, data.table(
      iteration = iter,
      a_male = a_m,
      b_male = b_m,
      a_female = a_f, 
      b_female = b_f,
      objective = objective
    ))
    
    # ========================================================================
    # Check convergence
    # ========================================================================
    param_change <- sqrt((a_m - a_m_old)^2 + (b_m - b_m_old)^2 + 
                        (a_f - a_f_old)^2 + (b_f - b_f_old)^2)
    
    if(param_change < tolerance) {
      cat("Converged after", iter, "iterations\n")
      return(list(
        a_male = a_m, b_male = b_m,
        a_female = a_f, b_female = b_f,
        converged = TRUE,
        iterations = iter,
        path = path
      ))
    }
  }
  
  cat("Did not converge after", max_iter, "iterations\n")
  return(list(
    a_male = a_m, b_male = b_m,
    a_female = a_f, b_female = b_f,
    converged = FALSE,
    iterations = max_iter,
    path = path
  ))
}

# ============================================================================
# Regularization Path Analysis
# ============================================================================
regularization_path <- function(dt, 
                               lambda_sequence = 10^seq(-3, 2, by = 0.5),
                               equal_penalties = TRUE) {
  
  cat("Computing regularization path with", length(lambda_sequence), "lambda values...\n")
  
  path_results <- data.table()
  
  for(i in seq_along(lambda_sequence)) {
    lambda <- lambda_sequence[i]
    
    # Use equal penalties for intercept and slope, or vary independently
    lambda_a <- lambda
    lambda_b <- if(equal_penalties) lambda else lambda * 0.1  # slope penalty could be different
    
    # Fit model
    fit <- fit_hierarchical_model(dt, lambda_a = lambda_a, lambda_b = lambda_b)
    
    # Store results
    path_results <- rbind(path_results, data.table(
      lambda = lambda,
      lambda_a = lambda_a,
      lambda_b = lambda_b,
      a_male = fit$a_male,
      b_male = fit$b_male,
      a_female = fit$a_female,
      b_female = fit$b_female,
      a_global = fit$a_global,
      b_global = fit$b_global,
      rmse_overall = fit$rmse_overall,
      r2_male = fit$r2_male,
      r2_female = fit$r2_female,
      # Parameter differences
      a_diff = abs(fit$a_male - fit$a_female),
      b_diff = abs(fit$b_male - fit$b_female),
      # Shrinkage toward global
      a_male_shrinkage = abs(fit$a_male - fit$a_global),
      a_female_shrinkage = abs(fit$a_female - fit$a_global),
      b_male_shrinkage = abs(fit$b_male - fit$b_global),
      b_female_shrinkage = abs(fit$b_female - fit$b_global)
    ))
    
    if(i %% 5 == 0) cat("  Completed", i, "of", length(lambda_sequence), "fits\n")
  }
  
  return(path_results)
}

# ============================================================================
# Cross-Validation Framework  
# ============================================================================
cross_validate <- function(dt, lambda_sequence = 10^seq(-2, 1, by = 0.5), k_folds = 5) {
  
  # Create folds stratified by gender
  set.seed(123)
  dt[, fold := sample(rep(1:k_folds, length.out = .N)), by = gender]
  
  cv_results <- data.table()
  
  cat("Running", k_folds, "-fold cross-validation...\n")
  
  for(lambda in lambda_sequence) {
    fold_errors <- numeric(k_folds)
    
    for(fold in 1:k_folds) {
      # Split data
      dt_train <- dt[fold != !!fold]
      dt_test <- dt[fold == !!fold]
      
      # Fit on training data
      fit <- fit_hierarchical_model(dt_train, lambda_a = lambda, lambda_b = lambda)
      
      # Predict on test data
      dt_test_m <- dt_test[gender == "Male"]
      dt_test_f <- dt_test[gender == "Female"]
      
      if(nrow(dt_test_m) > 0) {
        pred_m <- fit$a_male + fit$b_male * dt_test_m$weight_075
        error_m <- (dt_test_m$bmr - pred_m)^2
      } else {
        error_m <- numeric(0)
      }
      
      if(nrow(dt_test_f) > 0) {
        pred_f <- fit$a_female + fit$b_female * dt_test_f$weight_075  
        error_f <- (dt_test_f$bmr - pred_f)^2
      } else {
        error_f <- numeric(0)
      }
      
      fold_errors[fold] <- mean(c(error_m, error_f))
    }
    
    cv_results <- rbind(cv_results, data.table(
      lambda = lambda,
      cv_mse = mean(fold_errors),
      cv_rmse = sqrt(mean(fold_errors)),
      cv_se = sd(fold_errors) / sqrt(k_folds)
    ))
  }
  
  # Find optimal lambda
  optimal_idx <- which.min(cv_results$cv_mse)
  optimal_lambda <- cv_results$lambda[optimal_idx]
  
  cat("Optimal lambda:", optimal_lambda, "(CV RMSE:", 
      round(cv_results$cv_rmse[optimal_idx], 2), ")\n")
  
  return(list(
    cv_results = cv_results,
    optimal_lambda = optimal_lambda,
    data_with_folds = dt
  ))
}

# ============================================================================
# Visualization Functions
# ============================================================================
plot_results <- function(fit_result = NULL, path_result = NULL, cv_result = NULL) {
  
  plots <- list()
  
  # Plot 1: Parameter shrinkage path
  if(!is.null(path_result)) {
    path_long <- melt(path_result, 
                     id.vars = c("lambda"),
                     measure.vars = c("a_male", "a_female", "b_male", "b_female"),
                     variable.name = "parameter", value.name = "value")
    
    plots$shrinkage_path <- ggplot(path_long, aes(x = log10(lambda), y = value, 
                                                 color = parameter, linetype = parameter)) +
      geom_line(size = 1) +
      geom_hline(data = path_result[1], aes(yintercept = a_global), 
                linetype = "dashed", alpha = 0.7, color = "black") +
      geom_hline(data = path_result[1], aes(yintercept = b_global), 
                linetype = "dashed", alpha = 0.7, color = "black") +
      labs(title = "Regularization Path: Parameter Shrinkage",
           subtitle = "Dashed lines show global parameter values",
           x = "log10(λ)", y = "Parameter Value") +
      theme_minimal() +
      theme(legend.position = "bottom")
  }
  
  # Plot 2: Model performance vs regularization
  if(!is.null(path_result)) {
    plots$performance_path <- ggplot(path_result, aes(x = log10(lambda))) +
      geom_line(aes(y = rmse_overall, color = "RMSE"), size = 1) +
      geom_line(aes(y = (1 - r2_male) * 100, color = "1 - R² (Male)"), size = 1) +
      geom_line(aes(y = (1 - r2_female) * 100, color = "1 - R² (Female)"), size = 1) +
      labs(title = "Model Performance vs Regularization",
           x = "log10(λ)", y = "Performance Metric",
           color = "Metric") +
      theme_minimal() +
      theme(legend.position = "bottom")
  }
  
  # Plot 3: Cross-validation curve
  if(!is.null(cv_result)) {
    plots$cv_curve <- ggplot(cv_result$cv_results, aes(x = log10(lambda))) +
      geom_line(aes(y = cv_rmse), size = 1, color = "blue") +
      geom_ribbon(aes(ymin = cv_rmse - cv_se, ymax = cv_rmse + cv_se), 
                 alpha = 0.3, fill = "blue") +
      geom_vline(xintercept = log10(cv_result$optimal_lambda), 
                linetype = "dashed", color = "red") +
      labs(title = "Cross-Validation Curve",
           subtitle = paste("Optimal λ =", round(cv_result$optimal_lambda, 4)),
           x = "log10(λ)", y = "CV RMSE") +
      theme_minimal()
  }
  
  # Plot 4: Parameter differences
  if(!is.null(path_result)) {
    plots$parameter_differences <- ggplot(path_result, aes(x = log10(lambda))) +
      geom_line(aes(y = a_diff, color = "Intercept difference"), size = 1) +
      geom_line(aes(y = b_diff, color = "Slope difference"), size = 1) +
      labs(title = "Parameter Differences Between Genders",
           x = "log10(λ)", y = "Absolute Difference",
           color = "Parameter") +
      theme_minimal() +
      theme(legend.position = "bottom")
  }
  
  return(plots)
}

# ============================================================================
# Complete Analysis Runner
# ============================================================================
run_full_analysis <- function(n_points = 30, 
                             lambda_sequence = 10^seq(-2, 1.5, by = 0.3),
                             k_folds = 5) {
  
  cat("=== Running Complete Regularized BMR Analysis ===\n\n")
  
  # Generate data
  cat("1. Generating synthetic Harris-Benedict data...\n")
  dt <- generate_hb_data(n_points = n_points)
  
  # Fit single model
  cat("\n2. Fitting hierarchical model with default regularization...\n")
  fit <- fit_hierarchical_model(dt, lambda_a = 0.1, lambda_b = 0.1)
  
  # Regularization path
  cat("\n3. Computing regularization path...\n")
  path <- regularization_path(dt, lambda_sequence = lambda_sequence)
  
  # Cross-validation
  cat("\n4. Cross-validation to find optimal lambda...\n")
  cv <- cross_validate(dt, lambda_sequence = lambda_sequence, k_folds = k_folds)
  
  # Generate plots
  cat("\n5. Creating visualizations...\n")
  plots <- plot_results(fit_result = fit, path_result = path, cv_result = cv)
  
  # Summary
  cat("\n=== Analysis Summary ===\n")
  cat("Data points:", nrow(dt), "\n")
  cat("Optimal lambda:", cv$optimal_lambda, "\n")
  cat("Parameter ranges in path:\n")
  cat("  Intercept difference: [", round(min(path$a_diff), 2), ",", 
      round(max(path$a_diff), 2), "]\n")
  cat("  Slope difference: [", round(min(path$b_diff), 2), ",", 
      round(max(path$b_diff), 2), "]\n")
  
  # Return comprehensive results
  results <- list(
    data = dt,
    single_fit = fit,
    regularization_path = path,
    cross_validation = cv,
    plots = plots,
    
    # Summary statistics
    summary = list(
      n_observations = nrow(dt),
      optimal_lambda = cv$optimal_lambda,
      parameter_ranges = list(
        intercept_diff_range = range(path$a_diff),
        slope_diff_range = range(path$b_diff)
      )
    )
  )
  
  class(results) <- "bmr_full_analysis"
  return(results)
}

# ============================================================================
# Usage Examples and Testing
# ============================================================================

# Initialize the analysis framework
h <- bmr_regularized_analysis()
str(h)
# Quick analysis
results <- h$run_full_analysis()

# Access individual plots
results$plots$shrinkage_path
results$plots$cv_curve

# Custom analysis
dt <- h$generate_hb_data(n_points = 40)
fit <- h$fit_hierarchical_model(dt, lambda_a = 0.5, lambda_b = 0.1) 
path <- h$regularization_path(dt)

cat("BMR Regularized Regression Framework Loaded\n")
cat("Initialize with: h <- bmr_regularized_analysis()\n")
cat("Run analysis with: results <- h$run_full_analysis()\n")
