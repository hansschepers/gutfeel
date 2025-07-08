server <- function(input, output, session) {
  
  # Reactive data generation
  bmr_data <- reactive({
    # Get ranges from editable table if available
    bmi_range <- c(18, 35)
    height_range <- c(150, 180)
    age_range <- c(18, 70)
    
    if (!is.null(input$ranges_table)) {
      ranges_df <- hot_to_r(input$ranges_table)
      if (nrow(ranges_df) >= 3) {
        height_range <- c(ranges_df$Min_Value[2], ranges_df$Max_Value[2])
        age_range <- c(ranges_df$Min_Value[3], ranges_df$Max_Value[3])
        if (input$data_method == "bmi_driven" && nrow(ranges_df) >= 4) {
          bmi_range <- c(ranges_df$Min_Value[4], ranges_df$Max_Value[4])
        }
      }
    }
    
    # Get true parameters from combined table if available, otherwise use Harris-Benedict defaults
    true_betas <- HARRIS_BENEDICT$female  # Default to female since that's selected
    true_alphas <- c(1, 1, 1)
    
    # Try to get values from combined parameter table if it exists
    if (!is.null(input$combined_params_table)) {
      params_df <- hot_to_r(input$combined_params_table)
      if (nrow(params_df) == 7) {  # 4 betas + 3 alphas
        true_betas <- list(
          intercept = params_df$True_Value[1],
          weight = params_df$True_Value[2], 
          height = params_df$True_Value[3],
          age = params_df$True_Value[4]
        )
        true_alphas <- c(params_df$True_Value[5], params_df$True_Value[6], params_df$True_Value[7])
      }
    }
    
    generate_bmr_data(
      n_samples = if(is.null(input$sample_size)) 100 else input$sample_size,
      noise_level = if(is.null(input$noise_level)) 5 else input$noise_level / 100,
      data_method = input$data_method,
      bmi_range = bmi_range,
      height_range = height_range,
      age_range = age_range,
      true_betas = true_betas,
      true_alphas = true_alphas
    )
  })
  
  # Reactive model fitting
  model_results <- reactive({
    # Get exponents from combined parameter table
    exp_w <- 1
    exp_h <- 1 
    exp_a <- 1
    
    if (!is.null(input$combined_params_table)) {
      params_df <- hot_to_r(input$combined_params_table)
      if (nrow(params_df) == 7) {  # 4 betas + 3 alphas
        exp_w <- params_df$Fitting_Value[5]  # α₁
        exp_h <- params_df$Fitting_Value[6]  # α₂
        exp_a <- params_df$Fitting_Value[7]  # α₃
      }
    }
    
    exponents <- c(exp_w, exp_h, exp_a)
    
    # Handle checkboxGroupInput for centering correctly
    centering <- c(
      "weight" %in% input$centering_vars,
      "height" %in% input$centering_vars,
      "age" %in% input$centering_vars
    )
    
    fit_bmr_model(
      dt = bmr_data(),
      exponents = exponents,
      centering = centering,
      gender_subset = input$gender_subset
    )
  })
  
  # Data method display (compact, no empty lines)
  output$data_method_display <- renderText({
    method_text <- switch(input$data_method,
                          "independent" = "Independent weight & height (Unrealistic correlations)",
                          "bmi_driven" = paste0("BMI-driven generation: BMI & height → weight (Realistic correlations)", 
                                                if (!is.null(input$bmi_min) && !is.null(input$bmi_max)) 
                                                  paste0("\nBMI range: ", input$bmi_min, "-", input$bmi_max) else ""),
                          "user_specified" = "User-specified true model: Custom β and α values"
    )
    
    # Add green/black line legend
    paste0(method_text, "\n\nGreen = True model\nBlack = Fitted model")
  })
  
  # Fixed version of combined_params_table in server.R
  output$combined_params_table <- renderRHandsontable({
    # Get current fitted results
    results <- model_results()
    
    # Use Harris-Benedict defaults based on selected gender
    default_true_betas <- if (input$gender_subset == "male") {
      HARRIS_BENEDICT$male
    } else if (input$gender_subset == "female") {
      HARRIS_BENEDICT$female
    } else {
      # For "both", use female as default
      HARRIS_BENEDICT$female
    }
    
    # Default true values (fallback)
    default_true_values <- c(default_true_betas[["intercept"]], default_true_betas[["weight"]], 
                             default_true_betas[["height"]], default_true_betas[["age"]],
                             1.0, 1.0, 1.0)
    
    # Check if user has already edited the table and preserve those values
    if (!is.null(input$combined_params_table)) {
      existing_df <- hot_to_r(input$combined_params_table)
      if (nrow(existing_df) == 7) {
        # Use existing user-edited True_Value column
        true_values_to_use <- existing_df$True_Value
      } else {
        # Fallback to defaults if table structure is wrong
        true_values_to_use <- default_true_values
      }
    } else {
      # First time rendering - use defaults
      true_values_to_use <- default_true_values
    }
    
    # Also preserve existing Fitting_Value column (for alpha exponents)
    if (!is.null(input$combined_params_table)) {
      existing_df <- hot_to_r(input$combined_params_table)
      if (nrow(existing_df) == 7) {
        fitting_values_to_use <- existing_df$Fitting_Value
      } else {
        fitting_values_to_use <- c(NA, NA, NA, NA, 1.0, 1.0, 1.0)
      }
    } else {
      fitting_values_to_use <- c(NA, NA, NA, NA, 1.0, 1.0, 1.0)
    }
    
    df <- data.frame(
      Parameter = c("β₀ (Intercept)", "β₁ (Weight)", "β₂ (Height)", "β₃ (Age)", 
                    "α₁ (Weight exp)", "α₂ (Height exp)", "α₃ (Age exp)"),
      Type = c("Beta", "Beta", "Beta", "Beta", "Alpha", "Alpha", "Alpha"),
      True_Value = true_values_to_use,  # Use preserved/user-edited values
      Fitting_Value = fitting_values_to_use,  # Preserve alpha exponents too
      Fitted_Estimate = c(results$coefficients, rep(NA, 3)),  # Only betas get fitted
      Fitted_SE = c(results$standard_errors, rep(NA, 3)),    # Only betas have SEs
      stringsAsFactors = FALSE
    )
    
    # Make conditional columns read-only based on data method
    rht <- rhandsontable(df, height = 210, width = 600) |>
      hot_col("Parameter", readOnly = TRUE) |>
      hot_col("Type", readOnly = TRUE) |>
      hot_col("Fitted_Estimate", readOnly = TRUE) |>
      hot_col("Fitted_SE", readOnly = TRUE) |>
      hot_col("True_Value", type = "numeric", format = "0.000") |>
      hot_col("Fitting_Value", type = "numeric", format = "0.000")
    
    # Make True_Value read-only for non-user-specified methods
    if (input$data_method != "user_specified") {
      rht <- rht |> hot_col("True_Value", readOnly = TRUE)
    }
    
    rht
  })
  
  
  # # Combined parameter table (editable)
  # output$combined_params_table <- renderRHandsontable({
  #   # Get current fitted results
  #   results <- model_results()
  #   
  #   # Use Harris-Benedict defaults based on selected gender
  #   default_true_betas <- if (input$gender_subset == "male") {
  #     HARRIS_BENEDICT$male
  #   } else if (input$gender_subset == "female") {
  #     HARRIS_BENEDICT$female
  #   } else {
  #     # For "both", use female as default
  #     HARRIS_BENEDICT$female
  #   }
  #   
  #   df <- data.frame(
  #     Parameter = c("β₀ (Intercept)", "β₁ (Weight)", "β₂ (Height)", "β₃ (Age)", 
  #                   "α₁ (Weight exp)", "α₂ (Height exp)", "α₃ (Age exp)"),
  #     Type = c("Beta", "Beta", "Beta", "Beta", "Alpha", "Alpha", "Alpha"),
  #     True_Value = c(default_true_betas[["intercept"]], default_true_betas[["weight"]], 
  #                    default_true_betas[["height"]], default_true_betas[["age"]],
  #                    1.0, 1.0, 1.0),
  #     Fitting_Value = c(NA, NA, NA, NA, 1.0, 1.0, 1.0),  # Betas fitted, alphas user-set
  #     Fitted_Estimate = c(results$coefficients, rep(NA, 3)),  # Only betas get fitted
  #     Fitted_SE = c(results$standard_errors, rep(NA, 3)),    # Only betas have SEs
  #     stringsAsFactors = FALSE
  #   )
  #   
  #   # Make conditional columns read-only based on data method
  #   rht <- rhandsontable(df, height = 210, width = 600) |>
  #     hot_col("Parameter", readOnly = TRUE) |>
  #     hot_col("Type", readOnly = TRUE) |>
  #     hot_col("Fitted_Estimate", readOnly = TRUE) |>
  #     hot_col("Fitted_SE", readOnly = TRUE) |>
  #     hot_col("True_Value", type = "numeric", format = "0.000") |>
  #     hot_col("Fitting_Value", type = "numeric", format = "0.000")
  #   
  #   # Make True_Value read-only for non-user-specified methods
  #   if (input$data_method != "user_specified") {
  #     rht <- rht |> hot_col("True_Value", readOnly = TRUE)
  #   }
  #   
  #   rht
  # })
  
  # Editable ranges table (rhandsontable)
  output$ranges_table <- renderRHandsontable({
    # Current data to show actual ranges
    dt <- bmr_data()
    
    df <- data.frame(
      Predictor = c("Weight (kg)", "Height (cm)", "Age (years)", "BMI"),
      Min_Value = c(50, 150, 18, 18),
      Max_Value = c(90, 180, 70, 35),
      Actual_Range = c(
        sprintf("%.1f - %.1f", min(dt$weight), max(dt$weight)),
        sprintf("%.1f - %.1f", min(dt$height), max(dt$height)),
        sprintf("%.1f - %.1f", min(dt$age), max(dt$age)),
        if ("bmi" %in% names(dt)) sprintf("%.1f - %.1f", min(dt$bmi), max(dt$bmi)) else "Derived"
      ),
      stringsAsFactors = FALSE
    )
    
    # Make weight row read-only for BMI-driven, BMI row read-only for independent
    rht <- rhandsontable(df, height = 250, width = 450) |> 
      hot_col("Predictor", readOnly = TRUE) |>
      hot_col("Actual_Range", readOnly = TRUE) |>
      hot_col("Min_Value", type = "numeric", format = "0") |>
      hot_col("Max_Value", type = "numeric", format = "0")
    
    # Conditional read-only based on data method
    if (input$data_method == "bmi_driven") {
      rht <- rht |> hot_cell(1, "Min_Value", readOnly = TRUE) |>
        hot_cell(1, "Max_Value", readOnly = TRUE)  # Weight row
    } else if (input$data_method == "independent") {
      rht <- rht |> hot_cell(4, "Min_Value", readOnly = TRUE) |>
        hot_cell(4, "Max_Value", readOnly = TRUE)  # BMI row
    }
    
    rht
  })
  
  # Predictor correlation table (rhandsontable)
  output$predictor_correlation_table <- renderRHandsontable({
    dt <- bmr_data()
    pred_data <- dt[, .(weight, height, age)]
    cor_matrix <- cor(pred_data)
    
    # Convert correlation matrix to data frame
    cor_df <- data.frame(
      Variable = rownames(cor_matrix),
      Weight = round(cor_matrix[, "weight"], 3),
      Height = round(cor_matrix[, "height"], 3),
      Age = round(cor_matrix[, "age"], 3),
      stringsAsFactors = FALSE
    )
    
    rhandsontable(cor_df, height = 150, width = 350, readOnly = TRUE) |>
      hot_col("Variable", readOnly = TRUE)
  })
  
  # Parameter correlation table (rhandsontable)
  output$param_correlation_table <- renderRHandsontable({
    results <- model_results()
    param_names <- c("β₀", "β₁", "β₂", "β₃")
    
    cor_matrix <- results$correlation
    cor_df <- data.frame(
      Parameter = param_names,
      β0 = round(cor_matrix[, 1], 3),
      β1 = round(cor_matrix[, 2], 3),
      β2 = round(cor_matrix[, 3], 3),
      β3 = round(cor_matrix[, 4], 3),
      stringsAsFactors = FALSE
    )
    
    rhandsontable(cor_df, height = 220, width = 350, readOnly = TRUE) |>
      hot_col("Parameter", readOnly = TRUE)
  })
  
  # Parameter covariance table (rhandsontable)
  output$param_covariance_table <- renderRHandsontable({
    results <- model_results()
    param_names <- c("β₀", "β₁", "β₂", "β₃")
    
    cov_matrix <- results$covariance
    cov_df <- data.frame(
      Parameter = param_names,
      β0 = round(cov_matrix[, 1], 4),
      β1 = round(cov_matrix[, 2], 4),
      β2 = round(cov_matrix[, 3], 4),
      β3 = round(cov_matrix[, 4], 4),
      stringsAsFactors = FALSE
    )
    
    rhandsontable(cov_df, height = 220, width = 350, readOnly = TRUE) |>
      hot_col("Parameter", readOnly = TRUE)
  })
  
  # Predictor correlation plot
  output$predictor_corr_plot <- renderPlot({
    dt <- bmr_data()
    
    # Select only predictors for correlation
    pred_data <- dt[, .(weight, height, age)]
    cor_matrix <- cor(pred_data)
    
    corrplot(cor_matrix, method = "color", type = "full", order = "original",
             tl.cex = 0.9, tl.col = "black", tl.srt = 45,
             col = colorRampPalette(c("blue", "white", "red"))(100),
             addCoef.col = "black", number.cex = 0.8,
             title = "")
  })
  
  # Q-Q plot for residual normality
  output$qq_plot <- renderPlot({
    results <- model_results()
    
    if (is.null(results$residuals)) {
      return(ggplot() + geom_text(aes(x = 0, y = 0, label = "No residuals"), size = 5))
    }
    
    # Create Q-Q plot data
    residuals <- as.numeric(results$residuals)
    n <- length(residuals)
    theoretical_quantiles <- qnorm((1:n - 0.5) / n)
    sample_quantiles <- sort(residuals)
    
    qq_data <- data.frame(
      theoretical = theoretical_quantiles,
      sample = sample_quantiles
    )
    
    ggplot(qq_data, aes(x = theoretical, y = sample)) +
      geom_point(alpha = 0.6, size = 1.5) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      labs(x = "Theoretical Quantiles", y = "Sample Quantiles", 
           title = "Normal Q-Q Plot of Residuals") +
      theme_minimal() +
      theme(text = element_text(size = 10))
  })
  
  # Current model display
  output$current_model <- renderText({
    # Get exponents from combined table
    exp_w <- 1
    exp_h <- 1
    exp_a <- 1
    
    if (!is.null(input$combined_params_table)) {
      params_df <- hot_to_r(input$combined_params_table)
      if (nrow(params_df) == 7) {
        exp_w <- params_df$Fitting_Value[5]  # α₁
        exp_h <- params_df$Fitting_Value[6]  # α₂ 
        exp_a <- params_df$Fitting_Value[7]  # α₃
      }
    }
    
    sprintf("BMR = β₀ + β₁×weight^%.2f + β₂×height^%.2f + β₃×age^%.2f",
            exp_w, exp_h, exp_a)
  })
  
  # Model fit statistics
  output$fit_stats <- renderTable({
    results <- model_results()
    data.frame(
      Metric = c("R²", "MSE", "RMSE"),
      Value = c(
        sprintf("%.4f", results$r_squared),
        sprintf("%.2f", results$mse),
        sprintf("%.2f", sqrt(results$mse))
      )
    )
  }, colnames = FALSE, spacing = "xs", width = "100%")
  
  # Parameter estimates
  output$param_estimates <- renderTable({
    results <- model_results()
    param_names <- c("β₀", "β₁", "β₂", "β₃")
    
    # Check for large intercept with centering
    any_centered <- length(input$centering_vars) > 0
    intercept_large <- any_centered && abs(results$coefficients[1]) > 50
    
    data.frame(
      Parameter = param_names,
      Estimate = sprintf("%.3f", results$coefficients),
      SE = sprintf("%.3f", results$standard_errors),
      Warning = c(
        ifelse(intercept_large, "⚠️", ""),
        rep("", 3)
      )
    )
  }, colnames = TRUE, spacing = "xs", width = "100%")
  
  # Correlation plot
  output$correlation_plot <- renderPlot({
    results <- model_results()
    param_names <- c("β₀", "β₁", "β₂", "β₃")
    
    cor_mat <- results$correlation
    rownames(cor_mat) <- colnames(cor_mat) <- param_names
    
    corrplot(cor_mat, method = "color", type = "full", order = "original",
             tl.cex = 0.8, tl.col = "black", tl.srt = 45,
             col = colorRampPalette(c("blue", "white", "red"))(100),
             addCoef.col = "black", number.cex = 0.7)
  })
  
  # Covariance table  
  output$covariance_table <- DT::renderDataTable({
    results <- model_results()
    param_names <- c("β₀", "β₁", "β₂", "β₃")
    
    cov_dt <- data.table(results$covariance)
    setnames(cov_dt, param_names)
    cov_dt[, Parameter := param_names]
    setcolorder(cov_dt, "Parameter")
    
    # Round for display
    numeric_cols <- param_names
    cov_dt[, (numeric_cols) := lapply(.SD, function(x) round(x, 4)), .SDcols = numeric_cols]
    
    cov_dt
  }, options = list(pageLength = 4, dom = 't', scrollX = TRUE), rownames = FALSE)
  
  
  # Residual plot
  output$residual_plot <- renderPlot({
    results <- model_results()
    
    if (is.null(results$fitted) || is.null(results$residuals)) {
      return(ggplot() + geom_text(aes(x = 0, y = 0, label = "Model fitting error"), size = 5))
    }
    
    plot_dt <- data.frame(
      fitted = as.numeric(results$fitted),
      residuals = as.numeric(results$residuals)
    )
    
    ggplot(plot_dt, aes(x = fitted, y = residuals)) +
      geom_point(alpha = 0.6, size = 1) +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
      geom_smooth(method = "loess", se = FALSE, color = "blue", linewidth = 0.8) +
      labs(x = "Fitted BMR", y = "Residuals") +
      theme_minimal() +
      theme(text = element_text(size = 12))
  })
  
  
  # Helper function to get true model coefficients based on gender and method
  get_true_coeffs <- function(gender) {
    if (input$data_method == "user_specified") {
      return(list(
        intercept = if(is.null(input$true_intercept)) 1500 else input$true_intercept,
        weight = if(is.null(input$true_weight)) 15 else input$true_weight,
        height = if(is.null(input$true_height)) 5 else input$true_height,
        age = if(is.null(input$true_age)) -5 else input$true_age
      ))
    } else {
      if (gender == "male") {
        return(HARRIS_BENEDICT$male)
      } else if (gender == "female") {
        return(HARRIS_BENEDICT$female)
      } else {
        # For "both", use weighted average
        male_coeffs <- HARRIS_BENEDICT$male
        female_coeffs <- HARRIS_BENEDICT$female
        return(list(
          intercept = (male_coeffs[1] + female_coeffs[1]) / 2,
          weight = (male_coeffs[2] + female_coeffs[2]) / 2,
          height = (male_coeffs[3] + female_coeffs[3]) / 2,
          age = (male_coeffs[4] + female_coeffs[4]) / 2
        ))
      }
    }
  }
  
  # Helper functions for getting true coefficients and alphas
  get_true_coeffs_from_table <- function() {
    if (!is.null(input$combined_params_table)) {
      params_df <- hot_to_r(input$combined_params_table)
      if (nrow(params_df) == 7) {
        return(list(
          intercept = params_df$True_Value[1],
          weight = params_df$True_Value[2],
          height = params_df$True_Value[3],
          age = params_df$True_Value[4]
        ))
      }
    }
    
    # Fallback to Harris-Benedict based on gender
    if (input$gender_subset == "male") {
      return(HARRIS_BENEDICT$male)
    } else if (input$gender_subset == "female") {
      return(HARRIS_BENEDICT$female)
    } else {
      # For "both", use female as default
      return(HARRIS_BENEDICT$female)
    }
  }
  
  get_true_alphas_from_table <- function() {
    if (!is.null(input$combined_params_table)) {
      params_df <- hot_to_r(input$combined_params_table)
      if (nrow(params_df) == 7) {
        return(c(params_df$True_Value[5], params_df$True_Value[6], params_df$True_Value[7]))
      }
    }
    return(c(1, 1, 1))  # Default
  }
  
  # BMR vs Weight plot with larger fonts
  output$plot_weight <- renderPlot({
    results <- model_results()
    dt <- results$data
    
    # Get exponents from combined table
    exp_w <- 1
    exp_h <- 1
    exp_a <- 1
    
    if (!is.null(input$combined_params_table)) {
      params_df <- hot_to_r(input$combined_params_table)
      if (nrow(params_df) == 7) {
        exp_w <- params_df$Fitting_Value[5]  # α₁
        exp_h <- params_df$Fitting_Value[6]  # α₂
        exp_a <- params_df$Fitting_Value[7]  # α₃
      }
    }
    
    # Create prediction line based on current fitted model
    weight_range <- seq(min(dt$weight), max(dt$weight), length.out = 50)
    pred_data <- data.frame(
      weight = weight_range,
      height = mean(dt$height),
      age = mean(dt$age)
    )
    
    # Apply same transformations as in fitted model
    alpha1 <- exp_w + ifelse(exp_w == 0 | exp_w == 1, 0.0001, 0)
    alpha2 <- exp_h + ifelse(exp_h == 0 | exp_h == 1, 0.0001, 0)
    alpha3 <- exp_a + ifelse(exp_a == 0 | exp_a == 1, 0.0001, 0)
    
    weight_term <- sign(pred_data$weight) * abs(pred_data$weight + 1e-6)^alpha1
    height_term <- sign(pred_data$height) * abs(pred_data$height + 1e-6)^alpha2
    age_term <- sign(pred_data$age) * abs(pred_data$age + 1e-6)^alpha3
    
    X_pred <- as.matrix(cbind(1, weight_term, height_term, age_term))
    pred_data$bmr_fitted <- as.numeric(X_pred %*% results$coefficients)
    
    # Create true model line
    true_coeffs <- get_true_coeffs_from_table()
    
    if (input$data_method == "user_specified") {
      # Use true alphas for true model
      true_alphas <- get_true_alphas_from_table()
      
      true_alpha1 <- true_alphas[1] + ifelse(true_alphas[1] == 0 | true_alphas[1] == 1, 0.0001, 0)
      true_alpha2 <- true_alphas[2] + ifelse(true_alphas[2] == 0 | true_alphas[2] == 1, 0.0001, 0)
      true_alpha3 <- true_alphas[3] + ifelse(true_alphas[3] == 0 | true_alphas[3] == 1, 0.0001, 0)
      
      pred_data$bmr_true <- true_coeffs$intercept + 
        true_coeffs$weight * (pred_data$weight^true_alpha1) +
        true_coeffs$height * (pred_data$height^true_alpha2) + 
        true_coeffs$age * (pred_data$age^true_alpha3)
    } else {
      # Harris-Benedict (linear)
      pred_data$bmr_true <- true_coeffs[["intercept"]] + 
        true_coeffs[["weight"]] * pred_data$weight +
        true_coeffs[["height"]] * pred_data$height + 
        true_coeffs[["age"]] * pred_data$age
    }
    
    ggplot() +
      geom_point(data = dt, aes(x = weight, y = bmr, color = gender), alpha = 0.6, size = 2) +
      geom_line(data = pred_data, aes(x = weight, y = bmr_fitted), color = "black", linewidth = 1.5) +
      geom_line(data = pred_data, aes(x = weight, y = bmr_true), color = "green", linewidth = 1.2, alpha = 0.8) +
      labs(x = "Weight (kg)", y = "BMR (kcal/day)") +
      theme_minimal() +
      theme(legend.position = "bottom", 
            text = element_text(size = 12),
            axis.title = element_text(size = 13),
            axis.text = element_text(size = 13),
            legend.text = element_text(size = 11))
  })
  
  output$plot_height <- renderPlot({
    results <- model_results()
    dt <- results$data
    
    # Get exponents from combined table
    exp_w <- 1
    exp_h <- 1
    exp_a <- 1
    
    if (!is.null(input$combined_params_table)) {
      params_df <- hot_to_r(input$combined_params_table)
      if (nrow(params_df) == 7) {
        exp_w <- params_df$Fitting_Value[5]  # α₁
        exp_h <- params_df$Fitting_Value[6]  # α₂
        exp_a <- params_df$Fitting_Value[7]  # α₃
      }
    }
    
    # Create prediction line based on current fitted model
    height_range <- seq(min(dt$height), max(dt$height), length.out = 50)
    pred_data <- data.frame(
      weight = mean(dt$weight),
      height = height_range,
      age = mean(dt$age)
    )
    
    # Apply same transformations as in fitted model
    alpha1 <- exp_w + ifelse(exp_w == 0 | exp_w == 1, 0.0001, 0)
    alpha2 <- exp_h + ifelse(exp_h == 0 | exp_h == 1, 0.0001, 0)
    alpha3 <- exp_a + ifelse(exp_a == 0 | exp_a == 1, 0.0001, 0)
    
    weight_term <- sign(pred_data$weight) * abs(pred_data$weight + 1e-6)^alpha1
    height_term <- sign(pred_data$height) * abs(pred_data$height + 1e-6)^alpha2
    age_term <- sign(pred_data$age) * abs(pred_data$age + 1e-6)^alpha3
    
    X_pred <- as.matrix(cbind(1, weight_term, height_term, age_term))
    pred_data$bmr_fitted <- as.numeric(X_pred %*% results$coefficients)
    
    # Create true model line  
    true_coeffs <- get_true_coeffs_from_table()
    
    if (input$data_method == "user_specified") {
      # Use true alphas for true model
      true_alphas <- get_true_alphas_from_table()
      
      true_alpha1 <- true_alphas[1] + ifelse(true_alphas[1] == 0 | true_alphas[1] == 1, 0.0001, 0)
      true_alpha2 <- true_alphas[2] + ifelse(true_alphas[2] == 0 | true_alphas[2] == 1, 0.0001, 0)
      true_alpha3 <- true_alphas[3] + ifelse(true_alphas[3] == 0 | true_alphas[3] == 1, 0.0001, 0)
      
      pred_data$bmr_true <- true_coeffs$intercept + 
        true_coeffs$weight * (pred_data$weight^true_alpha1) +
        true_coeffs$height * (pred_data$height^true_alpha2) + 
        true_coeffs$age * (pred_data$age^true_alpha3)
    } else {
      # Harris-Benedict (linear)
      pred_data$bmr_true <- true_coeffs[["intercept"]] + 
        true_coeffs[["weight"]] * pred_data$weight +
        true_coeffs[["height"]] * pred_data$height + 
        true_coeffs[["age"]] * pred_data$age
    }
    
    ggplot() +
      geom_point(data = dt, aes(x = height, y = bmr, color = gender), alpha = 0.6, size = 2) +
      geom_line(data = pred_data, aes(x = height, y = bmr_fitted), color = "black", linewidth = 1.5) +
      geom_line(data = pred_data, aes(x = height, y = bmr_true), color = "green", linewidth = 1.2, alpha = 0.8) +
      labs(x = "Height (cm)", y = "BMR (kcal/day)") +
      theme_minimal() +
      theme(legend.position = "bottom", 
            text = element_text(size = 12),
            axis.title = element_text(size = 13),
            axis.text = element_text(size = 13),
            legend.text = element_text(size = 11))
  })
  
  output$plot_age <- renderPlot({
    results <- model_results()
    dt <- results$data
    
    # Get exponents from combined table
    exp_w <- 1
    exp_h <- 1
    exp_a <- 1
    
    if (!is.null(input$combined_params_table)) {
      params_df <- hot_to_r(input$combined_params_table)
      if (nrow(params_df) == 7) {
        exp_w <- params_df$Fitting_Value[5]  # α₁
        exp_h <- params_df$Fitting_Value[6]  # α₂
        exp_a <- params_df$Fitting_Value[7]  # α₃
      }
    }
    
    # Create prediction line based on current fitted model
    age_range <- seq(min(dt$age), max(dt$age), length.out = 50)
    pred_data <- data.frame(
      weight = mean(dt$weight),
      height = mean(dt$height),
      age = age_range
    )
    
    # Apply same transformations as in fitted model
    alpha1 <- exp_w + ifelse(exp_w == 0 | exp_w == 1, 0.0001, 0)
    alpha2 <- exp_h + ifelse(exp_h == 0 | exp_h == 1, 0.0001, 0)
    alpha3 <- exp_a + ifelse(exp_a == 0 | exp_a == 1, 0.0001, 0)
    
    weight_term <- sign(pred_data$weight) * abs(pred_data$weight + 1e-6)^alpha1
    height_term <- sign(pred_data$height) * abs(pred_data$height + 1e-6)^alpha2
    age_term <- sign(pred_data$age) * abs(pred_data$age + 1e-6)^alpha3
    
    X_pred <- as.matrix(cbind(1, weight_term, height_term, age_term))
    pred_data$bmr_fitted <- as.numeric(X_pred %*% results$coefficients)
    
    # Create true model line
    true_coeffs <- get_true_coeffs_from_table()
    
    if (input$data_method == "user_specified") {
      # Use true alphas for true model
      true_alphas <- get_true_alphas_from_table()
      
      true_alpha1 <- true_alphas[1] + ifelse(true_alphas[1] == 0 | true_alphas[1] == 1, 0.0001, 0)
      true_alpha2 <- true_alphas[2] + ifelse(true_alphas[2] == 0 | true_alphas[2] == 1, 0.0001, 0)
      true_alpha3 <- true_alphas[3] + ifelse(true_alphas[3] == 0 | true_alphas[3] == 1, 0.0001, 0)
      
      pred_data$bmr_true <- true_coeffs$intercept + 
        true_coeffs$weight * (pred_data$weight^true_alpha1) +
        true_coeffs$height * (pred_data$height^true_alpha2) + 
        true_coeffs$age * (pred_data$age^true_alpha3)
    } else {
      # Harris-Benedict (linear)  
      pred_data$bmr_true <- true_coeffs[["intercept"]] + 
        true_coeffs[["weight"]] * pred_data$weight +
        true_coeffs[["height"]] * pred_data$height + 
        true_coeffs[["age"]] * pred_data$age
    }
    
    ggplot() +
      geom_point(data = dt, aes(x = age, y = bmr, color = gender), alpha = 0.6, size = 2) +
      geom_line(data = pred_data, aes(x = age, y = bmr_fitted), color = "black", linewidth = 1.5) +
      geom_line(data = pred_data, aes(x = age, y = bmr_true), color = "green", linewidth = 1.2, alpha = 0.8) +
      labs(x = "Age (years)", y = "BMR (kcal/day)") +
      theme_minimal() +
      theme(legend.position = "bottom", 
            text = element_text(size = 12),
            axis.title = element_text(size = 13),
            axis.text = element_text(size = 13),
            legend.text = element_text(size = 11))
  })
  
  # Model diagnostics
  output$intercept_analysis <- renderText({
    results <- model_results()
    any_centered <- length(input$centering_vars) > 0
    intercept <- results$coefficients[1]
    se_intercept <- results$standard_errors[1]
    
    analysis_text <- paste0("Data method: ", input$data_method, " | Gender: ", input$gender_subset)
    
    if (any_centered) {
      analysis_text <- paste0(analysis_text, " | Centered: ", paste(input$centering_vars, collapse = ", "))
    }
    
    # Show data generation ranges
    if (input$data_method == "bmi_driven" && !is.null(input$bmi_min) && !is.null(input$bmi_max)) {
      analysis_text <- paste0(analysis_text, "\nBMI range: ", input$bmi_min, "-", input$bmi_max)
    }
    
    if (input$data_method != "independent" && !is.null(input$height_min) && !is.null(input$height_max)) {
      analysis_text <- paste0(analysis_text, " | Height: ", input$height_min, "-", input$height_max, " cm")
    }
    
    analysis_text <- paste0(analysis_text, "\n\n")
    
    if (!any_centered) {
      analysis_text <- paste0(analysis_text, "Center predictors to assess intercept significance. ")
    } else {
      if (abs(intercept) > 2 * se_intercept) {
        if (abs(intercept) > 50) {
          analysis_text <- paste0(analysis_text, sprintf("β₀ = %.2f ± %.2f ⚠️ LARGE intercept! Suggests model inadequacy.", 
                                                         intercept, se_intercept))
        } else {
          analysis_text <- paste0(analysis_text, sprintf("β₀ = %.2f ± %.2f ✓ Significant but reasonable.", intercept, se_intercept))
        }
      } else {
        analysis_text <- paste0(analysis_text, sprintf("β₀ = %.2f ± %.2f ✓ Not significantly different from 0.", intercept, se_intercept))
      }
    }
    
    # Add correlation structure insight
    dt <- bmr_data()
    weight_height_cor <- cor(dt$weight, dt$height)
    analysis_text <- paste0(analysis_text, sprintf("\n\nWeight-Height correlation: %.3f", weight_height_cor))
    
    if (input$data_method == "bmi_driven") {
      analysis_text <- paste0(analysis_text, " (BMI-driven: realistic)")
    } else if (input$data_method == "independent") {
      analysis_text <- paste0(analysis_text, " (Independent: may be unrealistic)")
    }
    .analysis_text <<- analysis_text
    analysis_text
  })
}
