# Helper functions for matrix operations
multiply_matrices <- function(A, B) {
  A %*% B
}

invert_matrix_safe <- function(X) {
  tryCatch({
    solve(X)
  }, error = function(e) {
    # Use pseudoinverse if singular
    svd_result <- svd(X)
    svd_result$v %*% diag(1/svd_result$d) %*% t(svd_result$u)
  })
}

# Enhanced data generation function
generate_bmr_data <- function(n_samples = 100, noise_level = 0.05, 
                              data_method = "independent",
                              bmi_range = c(18, 35), height_range = c(150, 180), age_range = c(18, 70),
                              true_betas = list(intercept = 1500, weight = 15, height = 5, age = -5),
                              true_alphas = c(1, 1, 1)) {
  set.seed(12345)  # Reproducible data
  
  # Generate realistic physiological data
  dt <- data.table(
    id = 1:n_samples,
    gender = rep(c("male", "female"), each = n_samples/2)
  )
  
  if (data_method == "independent") {
    # Original approach: independent weight and height (unrealistic)
    dt[gender == "male", `:=`(
      weight = runif(.N, 60, 95),          # kg - males
      height = runif(.N, height_range[1] + 10, height_range[2] + 10),  # cm - males  
      age = runif(.N, age_range[1], age_range[2])
    )]
    
    dt[gender == "female", `:=`(
      weight = runif(.N, 45, 80),          # kg - females
      height = runif(.N, height_range[1], height_range[2]),     # cm - females
      age = runif(.N, age_range[1], age_range[2])
    )]
    
  } else if (data_method == "bmi_driven") {
    # BMI-driven approach: BMI and height independent, weight derived
    dt[gender == "male", `:=`(
      bmi = runif(.N, bmi_range[1] + 2, bmi_range[2] + 2),    # Males slightly higher BMI
      height = runif(.N, height_range[1] + 10, height_range[2] + 10),
      age = runif(.N, age_range[1], age_range[2])
    )]
    
    dt[gender == "female", `:=`(
      bmi = runif(.N, bmi_range[1], bmi_range[2]),
      height = runif(.N, height_range[1], height_range[2]),
      age = runif(.N, age_range[1], age_range[2])
    )]
    
    # Calculate weight from BMI and height
    dt[, weight := bmi * (height/100)^2]
    
  } else if (data_method == "user_specified") {
    # User-specified model with user-defined true betas
    dt[gender == "male", `:=`(
      weight = runif(.N, 60, 95),
      height = runif(.N, height_range[1] + 10, height_range[2] + 10),
      age = runif(.N, age_range[1], age_range[2])
    )]
    
    dt[gender == "female", `:=`(
      weight = runif(.N, 45, 80),
      height = runif(.N, height_range[1], height_range[2]),
      age = runif(.N, age_range[1], age_range[2])
    )]
  }
  
  # Generate BMR based on method and true parameters
  if (data_method == "user_specified") {
    # Use user-specified true betas and alphas
    alpha1 <- true_alphas[1] + ifelse(true_alphas[1] == 0 | true_alphas[1] == 1, 0.0001, 0)
    alpha2 <- true_alphas[2] + ifelse(true_alphas[2] == 0 | true_alphas[2] == 1, 0.0001, 0)
    alpha3 <- true_alphas[3] + ifelse(true_alphas[3] == 0 | true_alphas[3] == 1, 0.0001, 0)
    
    dt[, bmr := true_betas$intercept + 
         true_betas$weight * (weight^alpha1) +
         true_betas$height * (height^alpha2) + 
         true_betas$age * (age^alpha3)]
    
  } else {
    # Use original Harris-Benedict equations for other methods
    dt[gender == "male", bmr := 
         HARRIS_BENEDICT$male["intercept"] + 
         HARRIS_BENEDICT$male["weight"] * weight +
         HARRIS_BENEDICT$male["height"] * height + 
         HARRIS_BENEDICT$male["age"] * age]
    
    dt[gender == "female", bmr := 
         HARRIS_BENEDICT$female["intercept"] + 
         HARRIS_BENEDICT$female["weight"] * weight +
         HARRIS_BENEDICT$female["height"] * height + 
         HARRIS_BENEDICT$female["age"] * age]
  }
  
  # Store true BMR before adding noise
  dt[, bmr_true := bmr]
  
  # Add noise (proper percentage of mean BMR)
  dt[, bmr := bmr + rnorm(.N, 0, noise_level * mean(bmr))]
  
  return(dt)
}

# Model fitting function
fit_bmr_model <- function(dt, exponents = c(1,1,1), 
                          centering = c(FALSE, FALSE, FALSE), gender_subset = "both") {
  
  # Filter data if needed
  if (gender_subset != "both") {
    dt <- dt[gender == gender_subset]
  }
  
  # Apply centering
  dt_work <- copy(dt)
  predictors <- c("weight", "height", "age")
  
  for (i in seq_along(predictors)) {
    if (centering[i]) {
      pred <- predictors[i]
      dt_work[, (pred) := get(pred) - mean(get(pred))]
    }
  }
  
  # Prepare design matrix (always nonlinear with exponents)
  alpha1 <- exponents[1] + ifelse(exponents[1] == 0 | exponents[1] == 1, 0.0001, 0)
  alpha2 <- exponents[2] + ifelse(exponents[2] == 0 | exponents[2] == 1, 0.0001, 0)
  alpha3 <- exponents[3] + ifelse(exponents[3] == 0 | exponents[3] == 1, 0.0001, 0)
  
  weight_term <- sign(dt_work$weight) * abs(dt_work$weight + 1e-6)^alpha1
  height_term <- sign(dt_work$height) * abs(dt_work$height + 1e-6)^alpha2
  age_term <- sign(dt_work$age) * abs(dt_work$age + 1e-6)^alpha3
  
  X <- as.matrix(cbind(1, weight_term, height_term, age_term))
  y <- dt_work$bmr
  
  # Fit model: beta = (X'X)^-1 X'y
  XtX <- t(X) %*% X
  XtX_inv <- invert_matrix_safe(XtX)
  Xty <- t(X) %*% y
  beta <- XtX_inv %*% Xty
  
  # Calculate residuals and MSE
  y_pred <- as.numeric(X %*% beta)
  residuals <- as.numeric(y - y_pred)
  n <- length(y)
  p <- ncol(X)
  mse <- sum(residuals^2) / max(1, n - p)
  
  # Parameter covariance matrix
  cov_matrix <- mse * XtX_inv
  
  # Correlation matrix
  se <- sqrt(diag(cov_matrix))
  cor_matrix <- cov_matrix / outer(se, se)
  
  # R-squared
  ss_res <- sum(residuals^2)
  ss_tot <- sum((y - mean(y))^2)
  r_squared <- 1 - ss_res/ss_tot
  
  return(list(
    coefficients = as.numeric(beta),
    covariance = cov_matrix,
    correlation = cor_matrix,
    standard_errors = as.numeric(se),
    residuals = residuals,
    fitted = y_pred,
    r_squared = r_squared,
    mse = mse,
    data = dt_work,
    design_matrix = X
  ))
}

