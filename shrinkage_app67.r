library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(data.table)
library(quantreg)
library(e1071)  # for SVR

# Generate power law BMR data with controllable outliers
generate_bmr_data <- function(n = 200, noise_sd = 50, exponent = 0.9, 
                              outlier_pct = 5, outlier_amplitude = 2, seed_value = 42) {
  if (seed_value != 0) {
    set.seed(seed_value)  # Only set seed if not 0
  }
  
  # Weight range that creates good curvature
  weight <- runif(n, 40, 120)  # kg, wide range for power law effect
  
  # Power law: BMR = 100 + 30 * weight^exponent (kcal/day)
  # Variable exponent creates different curvature effects
  bmr_true <- 100 + 30 * (weight^exponent)
  
  # Add Gaussian noise
  bmr_observed <- bmr_true + rnorm(n, 0, noise_sd)
  
  # Add systematic outliers (preferentially at high weights to force anti-shrinkage)
  n_outliers <- round(n * outlier_pct / 100)
  if (n_outliers > 0) {
    # Bias outlier selection toward high weights (where we want to force overprediction)
    weight_probs <- pmax(0, weight - median(weight))  # Higher weight = higher prob
    weight_probs <- weight_probs / sum(weight_probs)
    
    outlier_idx <- sample(n, n_outliers, prob = weight_probs, replace = FALSE)
    
    # Add positive outliers (high BMR values) with amplitude control
    outlier_boost <- outlier_amplitude * noise_sd * abs(rnorm(n_outliers, 2, 0.5))
    bmr_observed[outlier_idx] <- bmr_observed[outlier_idx] + outlier_boost
  }
  
  # Ensure BMR > 50 (physiological minimum)
  bmr_observed <- pmax(bmr_observed, 50)
  
  data.table(
    weight = weight,         # kg
    bmr_true = bmr_true,     # kcal/day (true power law)
    bmr_obs = bmr_observed,  # kcal/day (observed with noise + outliers)
    is_outlier = 1:n %in% if(exists("outlier_idx")) outlier_idx else integer(0)
  )
}

# Custom expectile regression
expectile_regression <- function(x, y, tau = 0.5) {
  # Asymmetric least squares for expectiles
  loss_fun <- function(beta) {
    fitted <- cbind(1, x) %*% beta
    residuals <- y - fitted
    weights <- ifelse(residuals >= 0, tau, 1 - tau)
    sum(weights * residuals^2)
  }
  
  # Start with OLS estimates
  ols_fit <- lm(y ~ x)
  beta_start <- coef(ols_fit)
  
  # Optimize
  result <- optim(beta_start, loss_fun, method = "BFGS")
  
  list(
    coefficients = result$par,
    fitted_values = cbind(1, x) %*% result$par,
    tau = tau
  )
}

# Custom asymmetric loss regression
asymmetric_loss_regression <- function(x, y, alpha = 2) {
  loss_fun <- function(beta) {
    fitted <- cbind(1, x) %*% beta
    residuals <- y - fitted
    # alpha > 1: penalize underprediction more heavily
    loss <- ifelse(residuals > 0, alpha * residuals^2, residuals^2)
    sum(loss)
  }
  
  ols_fit <- lm(y ~ x)
  beta_start <- coef(ols_fit)
  
  result <- optim(beta_start, loss_fun, method = "BFGS")
  
  list(
    coefficients = result$par,
    fitted_values = cbind(1, x) %*% result$par,
    alpha = alpha
  )
}

# Calculate slope of predicted vs observed
calc_pred_obs_slope <- function(y_obs, y_pred) {
  if (length(unique(y_pred)) == 1) return(0)
  lm(y_pred ~ y_obs)$coefficients[2]
}

# Calculate R-squared
calc_r_squared <- function(y_obs, y_pred) {
  ss_res <- sum((y_obs - y_pred)^2)
  ss_tot <- sum((y_obs - mean(y_obs))^2)
  1 - (ss_res / ss_tot)
}

# Calculate RMSE  
calc_rmse <- function(y_obs, y_pred) {
  sqrt(mean((y_obs - y_pred)^2))
}

# UI
ui <- dashboardPage(
  dashboardHeader(title = "Shrinkage vs Anti-Shrinkage Regression Methods"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Main Analysis", tabName = "main", icon = icon("chart-line")),
      menuItem("Hyperparameter Effects", tabName = "hyper", icon = icon("sliders-h")),
      menuItem("Loss Functions", tabName = "loss", icon = icon("function")),
      menuItem("Data Parameters & Overview", tabName = "data", icon = icon("table"))
    ),
    
    br(),
    h4("Method Parameters", style = "margin-left: 15px;"),
    
    # Quantile regression tau
    sliderInput("quantile_tau", 
                "Quantile τ (0.5 = median, 0.95 = 95th percentile)",
                min = 0.1, max = 0.99, value = 0.8, step = 0.01),
    
    # Expectile regression tau
    sliderInput("expectile_tau",
                "Expectile τ (0.5 = mean, >0.5 = anti-shrinkage)",
                min = 0.1, max = 0.99, value = 0.7, step = 0.01),
    
    # Asymmetric loss alpha
    sliderInput("asymmetric_alpha",
                "Asymmetric Loss α (1 = symmetric, >1 = penalize underprediction)",
                min = 0.5, max = 10, value = 2, step = 0.1),
    
    # SVR epsilon
    sliderInput("svr_epsilon",
                "SVR ε-tube (larger = more tolerance for small errors)",
                min = 10, max = 200, value = 50, step = 10),
    
    # SVR cost
    sliderInput("svr_cost",
                "SVR Cost (higher = force fit to outliers outside tube)",
                min = 1, max = 500, value = 100, step = 25),
    
    br(),
    actionButton("regenerate_data", "Regenerate Data", 
                 icon = icon("refresh"), class = "btn-warning")
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .slope-positive { color: #28a745; font-weight: bold; }
        .slope-negative { color: #dc3545; font-weight: bold; }
        .method-name { font-weight: bold; color: #007bff; }
      "))
    ),
    
    tabItems(
      # Main Analysis Tab
      tabItem(tabName = "main",
        fluidRow(
          box(title = "BMR vs Weight: Different Regression Methods", 
              status = "primary", solidHeader = TRUE, width = 12,
              plotOutput("main_plot", height = "500px")
          )
        ),
        
        fluidRow(
          box(title = "Predicted vs Observed (Slope Analysis)", 
              status = "info", solidHeader = TRUE, width = 8,
              plotOutput("pred_obs_plot", height = "450px")
          ),
          
          box(title = "Slope Summary", 
              status = "success", solidHeader = TRUE, width = 4,
              div(style = "font-size: 14px;",
                tableOutput("slope_table")
              ),
              br(),
              div(style = "background: #f8f9fa; padding: 10px; border-radius: 5px;",
                h5("Interpretation:", style = "color: #495057;"),
                p("• Slope < 1: Conservative (shrinkage)", style = "margin: 5px 0;"),
                p("• Slope ≈ 1: Perfect prediction", style = "margin: 5px 0;"), 
                p("• Slope > 1: Anti-shrinkage", style = "margin: 5px 0;"),
                p("• Higher slopes → better extreme value prediction", style = "margin: 5px 0;")
              )
          )
        )
      ),
      
      # Hyperparameter Effects Tab
      tabItem(tabName = "hyper",
        fluidRow(
          box(title = "Slope vs Quantile Tau", 
              status = "success", solidHeader = TRUE, width = 6,
              plotOutput("quantile_slope_plot", height = "300px"),
              p("Higher τ values fit upper quantiles → should increase slope", 
                style = "font-size: 12px; color: #666;")
          ),
          
          box(title = "Slope vs Expectile Tau", 
              status = "warning", solidHeader = TRUE, width = 6,
              plotOutput("expectile_slope_plot", height = "300px"),
              p("τ > 0.5 penalizes underprediction → should increase slope", 
                style = "font-size: 12px; color: #666;")
          )
        ),
        
        fluidRow(
          box(title = "Slope vs Asymmetric Alpha", 
              status = "danger", solidHeader = TRUE, width = 6,
              plotOutput("asymmetric_slope_plot", height = "300px"),
              p("Higher α heavily penalizes underprediction → should increase slope", 
                style = "font-size: 12px; color: #666;")
          ),
          
          box(title = "Slope vs SVR Epsilon", 
              status = "primary", solidHeader = TRUE, width = 6,
              plotOutput("svr_slope_plot", height = "300px"),
              p("Larger ε ignores small errors → can lead to different behavior", 
                style = "font-size: 12px; color: #666;")
          )
        ),
        
        fluidRow(
          box(title = "How to Achieve Slope > 1?", 
              status = "info", solidHeader = TRUE, width = 12,
              div(style = "font-size: 14px; line-height: 1.6;",
                h4("Strategies for Anti-Shrinkage (Slope > 1):", style = "color: #17a2b8;"),
                tags$ul(
                  tags$li("Use extreme quantiles: τ = 0.95 or 0.99"),
                  tags$li("Use high asymmetric loss: α = 5-10"), 
                  tags$li("Use subcritical exponents: 0.7-0.9 (concave curves)"),
                  tags$li("Use moderate sample sizes: n = 100-300 (balance noise vs signal)"),
                  tags$li("Add strategic outliers: 5-15% with amplitude 2-4"),
                  tags$li("SVR: high cost C = 200-500 with moderate ε = 50-100")
                ),
                div(style = "background: #d1ecf1; padding: 10px; border-radius: 5px; margin-top: 10px;",
                    strong("Key Insight: "), 
                    "Subcritical exponents (< 1) create concave curves where linear fits ",
                    "systematically underestimate high values, enabling anti-shrinkage methods ",
                    "to achieve slope > 1 by chasing strategic outliers."
                )
              )
          )
        ),
        
        fluidRow(
          box(title = "How to Achieve Slope > 1?", 
              status = "info", solidHeader = TRUE, width = 12,
              div(style = "font-size: 14px; line-height: 1.6;",
                h4("Strategies for Anti-Shrinkage (Slope > 1):", style = "color: #17a2b8;"),
                tags$ul(
                  tags$li("Use extreme quantiles: τ = 0.95 or 0.99"),
                  tags$li("Use high asymmetric loss: α = 5-10"), 
                  tags$li("Use subcritical exponents: 0.7-0.9 (concave curves)"),
                  tags$li("Use moderate sample sizes: n = 100-300 (balance noise vs signal)"),
                  tags$li("Add strategic outliers: 5-15% with amplitude 2-4"),
                  tags$li("SVR: high cost C = 200-500 with moderate ε = 50-100")
                ),
                div(style = "background: #d1ecf1; padding: 10px; border-radius: 5px; margin-top: 10px;",
                    strong("Key Insight: "), 
                    "Subcritical exponents (< 1) create concave curves where linear fits ",
                    "systematically underestimate high values, enabling anti-shrinkage methods ",
                    "to achieve slope > 1 by chasing strategic outliers."
                )
              )
          )
        )
      ),
      
      # Loss Functions Tab
      tabItem(tabName = "loss",
        fluidRow(
          box(title = "Loss Function Equations", 
              status = "warning", solidHeader = TRUE, width = 12,
              div(style = "font-size: 16px; line-height: 1.8;",
                
                div(style = "background: #e3f2fd; padding: 15px; margin: 10px 0; border-radius: 5px;",
                    h4("True Model: BMR = 100 + 30 × weight^exponent", style = "color: #1565c0; margin: 0;"),
                    p("We fit STRAIGHT LINES to this curved relationship - this creates shrinkage!", style = "margin: 5px 0 0 0;"),
                    p("Exponent ≠ 1 creates curvature. Try exponent = 0.5 or 1.5 for dramatic effects!", style = "margin: 5px 0 0 0; font-style: italic;")
                ),
                
                h4("1. Ordinary Least Squares (OLS)", style = "color: #007bff;"),
                div(style = "background: #f8f9fa; padding: 15px; margin: 10px 0; border-left: 4px solid #007bff;",
                    p("L(y, ŷ) = (y - ŷ)²"),
                    p("Symmetric loss → maximum shrinkage toward mean")
                ),
                
                h4("2. Quantile Regression", style = "color: #28a745;"),
                div(style = "background: #f8f9fa; padding: 15px; margin: 10px 0; border-left: 4px solid #28a745;",
                    p("L(y, ŷ) = (y - ŷ) × [τ - I(y < ŷ)]"),
                    p("τ = quantile level (e.g., 0.8 = 80th percentile)"),
                    p("τ > 0.5 → fits upper envelope, reduces shrinkage for high values")
                ),
                
                h4("3. Expectile Regression", style = "color: #ffc107;"),
                div(style = "background: #f8f9fa; padding: 15px; margin: 10px 0; border-left: 4px solid #ffc107;",
                    p("L(y, ŷ) = |τ - I(y < ŷ)| × (y - ŷ)²"),
                    p("Asymmetric squared loss"),
                    p("τ > 0.5 → penalizes underprediction more → anti-shrinkage")
                ),
                
                h4("4. Asymmetric Loss", style = "color: #dc3545;"),
                div(style = "background: #f8f9fa; padding: 15px; margin: 10px 0; border-left: 4px solid #dc3545;",
                    p("L(y, ŷ) = α × (y - ŷ)² if y > ŷ, else (y - ŷ)²"),
                    p("α > 1 → heavily penalize underprediction"),
                    p("Forces straight line toward upper bounds of curved data")
                ),
                
                h4("5. Support Vector Regression", style = "color: #6f42c1;"),
                div(style = "background: #f8f9fa; padding: 15px; margin: 10px 0; border-left: 4px solid #6f42c1;",
                    p("L(y, ŷ) = max(0, |y - ŷ| - ε) + C × penalty"),
                    p("ε-insensitive loss → ignores small errors inside tube"),
                    p("Cost C → heavily penalizes large errors outside tube"),
                    p("High C + moderate ε → forces line toward outliers → anti-shrinkage")
                )
              )
          )
        )
      ),
      
      # Data Parameters & Overview Tab
      tabItem(tabName = "data",
        fluidRow(
          box(title = "Data Generation Parameters", 
              status = "warning", solidHeader = TRUE, width = 4,
              
              # Sample size
              numericInput("n_points", 
                           "Number of Data Points",
                           value = 200, min = 20, max = 1000, step = 20),
              
              # Power law exponent
              numericInput("exponent", 
                           "Power Law Exponent (BMR ∝ weight^exp)",
                           value = 0.9, min = 0.5, max = 1.8, step = 0.05),
              
              # Noise control
              numericInput("noise_sd", 
                           "Noise Standard Deviation (kcal/day)",
                           value = 50, min = 10, max = 300, step = 10),
              
              # Seed control
              numericInput("seed_value",
                           "Random Seed (0 = no reset)",
                           value = 42, min = 0, max = 9999, step = 1),
              
              br(),
              actionButton("regenerate_data", "Regenerate Data", 
                           icon = icon("refresh"), class = "btn-warning")
          ),
          
          box(title = "Outlier Controls", 
              status = "danger", solidHeader = TRUE, width = 4,
              
              # Outlier percentage
              sliderInput("outlier_pct", 
                          "Outlier Percentage (%)",
                          min = 0, max = 20, value = 5, step = 1),
              
              # Outlier amplitude
              sliderInput("outlier_amplitude", 
                          "Outlier Amplitude (× noise)",
                          min = 1, max = 5, value = 2, step = 0.2)
          ),
          
          box(title = "How to Achieve Slope > 1", 
              status = "info", solidHeader = TRUE, width = 4,
              div(style = "font-size: 13px; line-height: 1.4;",
                h5("🎯 Optimal Settings:", style = "color: #17a2b8; margin-top: 0;"),
                tags$ul(style = "margin: 8px 0; padding-left: 18px;",
                  tags$li("n = 100-300 (balance signal/noise)"),
                  tags$li("exponent = 0.7-0.9 (concave curves)"),
                  tags$li("outliers = 8-15% with amplitude 3-4"),
                  tags$li("quantile τ = 0.95-0.99"),
                  tags$li("asymmetric α = 5-10"),
                  tags$li("SVR: cost = 200-500, ε = 50-100")
                ),
                div(style = "background: #d1ecf1; padding: 8px; border-radius: 4px; margin-top: 8px;",
                    p(style = "margin: 0; font-size: 12px;",
                      strong("Key: "), "Concave curves + strategic outliers + extreme parameters = anti-shrinkage")
                )
              )
          )
        ),
        
        fluidRow(
          box(title = "Generated Power Law Dataset", 
              status = "info", solidHeader = TRUE, width = 12,
              div(style = "margin-bottom: 15px;",
                h4("True Power Law Relationship:"),
                div(style = "background: #e3f2fd; padding: 10px; border-radius: 5px; font-family: monospace; font-size: 16px;",
                    textOutput("power_law_equation")
                ),
                br(),
                p("We then fit STRAIGHT LINES to this curved relationship using different loss functions."),
                p("The power law creates curvature that linear models must approximate, leading to shrinkage."),
                p("Anti-shrinkage methods try to better capture the extremes of this curved relationship."),
                p(strong("Outliers (shown in red) are added preferentially at high weights to force anti-shrinkage behavior."))
              ),
              DT::dataTableOutput("data_table")
          )
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive data
  bmr_data <- reactive({
    input$regenerate_data  # Trigger regeneration
    generate_bmr_data(n = input$n_points, 
                      noise_sd = input$noise_sd, 
                      exponent = input$exponent,
                      outlier_pct = input$outlier_pct,
                      outlier_amplitude = input$outlier_amplitude,
                      seed_value = input$seed_value)
  })
  
  # Also regenerate when parameters change
  observeEvent(c(input$noise_sd, input$n_points, input$exponent, 
                 input$outlier_pct, input$outlier_amplitude, input$seed_value, input$svr_cost), {
    # Regenerate with new parameters
  })
  
  # Power law equation display
  output$power_law_equation <- renderText({
    paste0("BMR (kcal/day) = 100 + 30 × weight^", input$exponent)
  })
  
  # Reactive calculations
  regression_results <- reactive({
    # Force dependency on all inputs
    tau_q <- input$quantile_tau
    tau_e <- input$expectile_tau  
    alpha_a <- input$asymmetric_alpha
    eps_svr <- input$svr_epsilon
    cost_svr <- input$svr_cost
    
    data <- bmr_data()
    x <- data$weight  # kg
    y <- data$bmr_obs # kcal/day
    
    # 1. OLS
    ols_fit <- lm(y ~ x)
    
    # 2. Quantile regression
    qr_fit <- rq(y ~ x, tau = tau_q)
    
    # 3. Expectile regression
    exp_fit <- expectile_regression(x, y, tau = tau_e)
    
    # 4. Asymmetric loss
    asym_fit <- asymmetric_loss_regression(x, y, alpha = alpha_a)
    
    # 5. SVR (enhanced with higher cost for anti-shrinkage)
    svr_fit <- svm(x, y, type = "eps-regression", epsilon = eps_svr, 
                   cost = cost_svr, kernel = "linear", scale = FALSE)
    
    # Predictions for plotting
    x_pred <- seq(min(x), max(x), length.out = 100)
    
    list(
      data = data,
      x = x,
      y = y,
      x_pred = x_pred,
      
      ols = list(
        fit = ols_fit,
        pred = predict(ols_fit, newdata = data.frame(x = x_pred)),
        fitted = fitted(ols_fit)
      ),
      
      quantile = list(
        fit = qr_fit,
        pred = predict(qr_fit, newdata = data.frame(x = x_pred)),
        fitted = fitted(qr_fit)
      ),
      
      expectile = list(
        fit = exp_fit,
        pred = cbind(1, x_pred) %*% exp_fit$coefficients,
        fitted = exp_fit$fitted_values
      ),
      
      asymmetric = list(
        fit = asym_fit,
        pred = cbind(1, x_pred) %*% asym_fit$coefficients,
        fitted = asym_fit$fitted_values
      ),
      
      svr = list(
        fit = svr_fit,
        pred = predict(svr_fit, newdata = x_pred),
        fitted = predict(svr_fit, newdata = x)
      )
    )
  })
  
  # Main plot
  output$main_plot <- renderPlot({
    results <- regression_results()
    
    # Create named color vector
    vvv <- c("#007bff", "#28a745", "#ffc107", "#dc3545", "#6f42c1")
    names(vvv) <- c("OLS",
                    paste0("Quantile (τ=", input$quantile_tau, ")"),
                    paste0("Expectile (τ=", input$expectile_tau, ")"), 
                    paste0("Asymmetric (α=", input$asymmetric_alpha, ")"),
                    paste0("SVR (ε=", input$svr_epsilon, ", C=", input$svr_cost, ")"))
    
    # Create data frames for each line
    pred_data <- data.frame(
      x = rep(results$x_pred, 5),
      y = c(results$ols$pred, results$quantile$pred, results$expectile$pred,
            results$asymmetric$pred, results$svr$pred),
      method = rep(names(vvv), each = length(results$x_pred))
    )
    
    # Add normal points
    normal_data <- results$data[is_outlier == FALSE]
    outlier_data <- results$data[is_outlier == TRUE]
    
    p <- ggplot() +
      geom_point(data = normal_data, aes(x = weight, y = bmr_obs),
                 alpha = 0.6, size = 2, color = "gray30") +
      geom_point(data = outlier_data, aes(x = weight, y = bmr_obs),
                 alpha = 0.8, size = 3, color = "red", shape = 17) +  # triangles for outliers
      
      geom_line(data = pred_data, aes(x = x, y = y, color = method), 
                linewidth = 1.2) +
      
      scale_color_manual(values = vvv, name = "Method") +
      
      labs(
        title = "BMR vs Weight: Comparing Shrinkage vs Anti-Shrinkage Methods",
        x = "Weight (kg)",
        y = "BMR (kcal/day)",
        color = "Method"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom",
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10))
    
    print(p)
  })
  
  # Predicted vs Observed plot
  output$pred_obs_plot <- renderPlot({
    results <- regression_results()
    
    # Calculate slopes
    slopes <- c(
      ols = calc_pred_obs_slope(results$y, results$ols$fitted),
      quantile = calc_pred_obs_slope(results$y, results$quantile$fitted),
      expectile = calc_pred_obs_slope(results$y, results$expectile$fitted),
      asymmetric = calc_pred_obs_slope(results$y, results$asymmetric$fitted),
      svr = calc_pred_obs_slope(results$y, results$svr$fitted)
    )
    
    # Create named color vector
    vvv <- c("#007bff", "#28a745", "#ffc107", "#dc3545", "#6f42c1")
    names(vvv) <- c("OLS",
                    paste0("Quantile (τ=", input$quantile_tau, ")"),
                    paste0("Expectile (τ=", input$expectile_tau, ")"), 
                    paste0("Asymmetric (α=", input$asymmetric_alpha, ")"),
                    paste0("SVR (ε=", input$svr_epsilon, ", C=", input$svr_cost, ")"))
    
    # Create data for plotting
    pred_obs_data <- data.frame(
      observed = rep(results$y, 5),
      predicted = c(results$ols$fitted, results$quantile$fitted, 
                   results$expectile$fitted, results$asymmetric$fitted, 
                   results$svr$fitted),
      method = rep(names(vvv), each = length(results$y)),
      slope = rep(round(slopes, 3), each = length(results$y))
    )
    
    p <- ggplot() +
      geom_point(data = pred_obs_data, aes(x = observed, y = predicted, color = method), 
                 alpha = 0.6) +
      geom_smooth(data = pred_obs_data, aes(x = observed, y = predicted, color = method), 
                  method = "lm", se = FALSE, linewidth = 1.2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
      
      scale_color_manual(values = vvv) +
      
      facet_wrap(~method, scales = "free") +
      labs(
        title = "Predicted vs Observed: Slope Analysis",
        subtitle = "Dashed line = perfect prediction (slope = 1)",
        x = "Observed BMR (kcal/day)",
        y = "Predicted BMR (kcal/day)"
      ) +
      theme_minimal() +
      theme(legend.position = "none")
    
    print(p)
  })
  
  # Slope summary table with R² and RMSE
  output$slope_table <- renderTable({
    results <- regression_results()
    
    slopes <- c(
      calc_pred_obs_slope(results$y, results$ols$fitted),
      calc_pred_obs_slope(results$y, results$quantile$fitted),
      calc_pred_obs_slope(results$y, results$expectile$fitted),
      calc_pred_obs_slope(results$y, results$asymmetric$fitted),
      calc_pred_obs_slope(results$y, results$svr$fitted)
    )
    
    r_squared <- c(
      calc_r_squared(results$y, results$ols$fitted),
      calc_r_squared(results$y, results$quantile$fitted),
      calc_r_squared(results$y, results$expectile$fitted),
      calc_r_squared(results$y, results$asymmetric$fitted),
      calc_r_squared(results$y, results$svr$fitted)
    )
    
    rmse <- c(
      calc_rmse(results$y, results$ols$fitted),
      calc_rmse(results$y, results$quantile$fitted),
      calc_rmse(results$y, results$expectile$fitted),
      calc_rmse(results$y, results$asymmetric$fitted),
      calc_rmse(results$y, results$svr$fitted)
    )
    
    method_names <- c("OLS",
                     paste0("Quantile (τ=", input$quantile_tau, ")"),
                     paste0("Expectile (τ=", input$expectile_tau, ")"),
                     paste0("Asymmetric (α=", input$asymmetric_alpha, ")"),
                     paste0("SVR (ε=", input$svr_epsilon, ", C=", input$svr_cost, ")"))
    
    slope_data <- data.frame(
      Method = method_names,
      Slope = round(slopes, 3),
      R2 = round(r_squared, 3),
      RMSE = round(rmse, 1),
      Behavior = ifelse(slopes < 0.95, "Conservative", 
                       ifelse(slopes > 1.05, "Anti-Shrinkage", "Balanced"))
    )
    
    slope_data
  }, rownames = FALSE, striped = TRUE)
  
  # Hyperparameter analysis plots
  
  # Quantile tau analysis with dual y-axis
  output$quantile_slope_plot <- renderPlot({
    data <- bmr_data()
    x <- data$weight
    y <- data$bmr_obs
    
    tau_values <- seq(0.1, 0.99, by = 0.02)
    slopes <- sapply(tau_values, function(tau) {
      qr_fit <- rq(y ~ x, tau = tau)
      fitted_vals <- fitted(qr_fit)
      calc_pred_obs_slope(y, fitted_vals)
    })
    
    r_squared <- sapply(tau_values, function(tau) {
      qr_fit <- rq(y ~ x, tau = tau)
      fitted_vals <- fitted(qr_fit)
      calc_r_squared(y, fitted_vals)
    })
    
    slope_data <- data.frame(tau = tau_values, slope = slopes, r_squared = r_squared)
    
    # Scale R² to fit nicely with slope
    r2_scaled <- slope_data$r_squared * (max(slopes) - min(slopes)) / (max(r_squared) - min(r_squared)) + min(slopes)
    
    ggplot(slope_data, aes(x = tau)) +
      geom_line(aes(y = slope, color = "Slope"), linewidth = 1.5) +
      geom_line(aes(y = r2_scaled, color = "R²"), linewidth = 1.5, linetype = "dashed") +
      geom_point(aes(y = slope, color = "Slope"), size = 1.5, alpha = 0.7) +
      geom_point(aes(y = r2_scaled, color = "R²"), size = 1.5, alpha = 0.7) +
      geom_hline(yintercept = 1, linetype = "dotted", color = "red", alpha = 0.7) +
      geom_vline(xintercept = input$quantile_tau, linetype = "dotted", color = "blue") +
      scale_color_manual(values = c("Slope" = "#28a745", "R²" = "#17a2b8")) +
      scale_y_continuous(
        name = "Slope (Pred vs Obs)",
        sec.axis = sec_axis(
          trans = ~ (. - min(slopes)) / (max(slopes) - min(slopes)) * (max(r_squared) - min(r_squared)) + min(r_squared),
          name = "R² (Goodness of Fit)"
        )
      ) +
      labs(title = "Quantile Regression: Slope & R² vs τ",
           x = "Quantile τ", color = "Metric") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  # Expectile tau analysis with dual y-axis
  output$expectile_slope_plot <- renderPlot({
    data <- bmr_data()
    x <- data$weight
    y <- data$bmr_obs
    
    tau_values <- seq(0.1, 0.99, by = 0.02)
    slopes <- sapply(tau_values, function(tau) {
      exp_fit <- expectile_regression(x, y, tau = tau)
      fitted_vals <- exp_fit$fitted_values
      calc_pred_obs_slope(y, fitted_vals)
    })
    
    r_squared <- sapply(tau_values, function(tau) {
      exp_fit <- expectile_regression(x, y, tau = tau)
      fitted_vals <- exp_fit$fitted_values
      calc_r_squared(y, fitted_vals)
    })
    
    slope_data <- data.frame(tau = tau_values, slope = slopes, r_squared = r_squared)
    
    # Scale R² to fit nicely with slope
    r2_scaled <- slope_data$r_squared * (max(slopes) - min(slopes)) / (max(r_squared) - min(r_squared)) + min(slopes)
    
    ggplot(slope_data, aes(x = tau)) +
      geom_line(aes(y = slope, color = "Slope"), linewidth = 1.5) +
      geom_line(aes(y = r2_scaled, color = "R²"), linewidth = 1.5, linetype = "dashed") +
      geom_point(aes(y = slope, color = "Slope"), size = 1.5, alpha = 0.7) +
      geom_point(aes(y = r2_scaled, color = "R²"), size = 1.5, alpha = 0.7) +
      geom_hline(yintercept = 1, linetype = "dotted", color = "red", alpha = 0.7) +
      geom_vline(xintercept = input$expectile_tau, linetype = "dotted", color = "blue") +
      scale_color_manual(values = c("Slope" = "#ffc107", "R²" = "#17a2b8")) +
      scale_y_continuous(
        name = "Slope (Pred vs Obs)",
        sec.axis = sec_axis(
          trans = ~ (. - min(slopes)) / (max(slopes) - min(slopes)) * (max(r_squared) - min(r_squared)) + min(r_squared),
          name = "R² (Goodness of Fit)"
        )
      ) +
      labs(title = "Expectile Regression: Slope & R² vs τ",
           x = "Expectile τ", color = "Metric") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  # Asymmetric alpha analysis with dual y-axis
  output$asymmetric_slope_plot <- renderPlot({
    data <- bmr_data()
    x <- data$weight
    y <- data$bmr_obs
    
    alpha_values <- seq(0.5, 10, by = 0.3)
    slopes <- sapply(alpha_values, function(alpha) {
      asym_fit <- asymmetric_loss_regression(x, y, alpha = alpha)
      fitted_vals <- asym_fit$fitted_values
      calc_pred_obs_slope(y, fitted_vals)
    })
    
    r_squared <- sapply(alpha_values, function(alpha) {
      asym_fit <- asymmetric_loss_regression(x, y, alpha = alpha)
      fitted_vals <- asym_fit$fitted_values
      calc_r_squared(y, fitted_vals)
    })
    
    slope_data <- data.frame(alpha = alpha_values, slope = slopes, r_squared = r_squared)
    
    # Scale R² to fit nicely with slope
    r2_scaled <- slope_data$r_squared * (max(slopes) - min(slopes)) / (max(r_squared) - min(r_squared)) + min(slopes)
    
    ggplot(slope_data, aes(x = alpha)) +
      geom_line(aes(y = slope, color = "Slope"), linewidth = 1.5) +
      geom_line(aes(y = r2_scaled, color = "R²"), linewidth = 1.5, linetype = "dashed") +
      geom_point(aes(y = slope, color = "Slope"), size = 1.5, alpha = 0.7) +
      geom_point(aes(y = r2_scaled, color = "R²"), size = 1.5, alpha = 0.7) +
      geom_hline(yintercept = 1, linetype = "dotted", color = "red", alpha = 0.7) +
      geom_vline(xintercept = input$asymmetric_alpha, linetype = "dotted", color = "blue") +
      scale_color_manual(values = c("Slope" = "#dc3545", "R²" = "#17a2b8")) +
      scale_y_continuous(
        name = "Slope (Pred vs Obs)",
        sec.axis = sec_axis(
          trans = ~ (. - min(slopes)) / (max(slopes) - min(slopes)) * (max(r_squared) - min(r_squared)) + min(r_squared),
          name = "R² (Goodness of Fit)"
        )
      ) +
      labs(title = "Asymmetric Loss: Slope & R² vs α",
           x = "Asymmetric Loss α", color = "Metric") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  # SVR epsilon analysis with dual y-axis  
  output$svr_slope_plot <- renderPlot({
    data <- bmr_data()
    x <- data$weight
    y <- data$bmr_obs
    
    eps_values <- seq(10, 200, by = 15)
    slopes <- sapply(eps_values, function(eps) {
      tryCatch({
        svr_fit <- svm(x, y, type = "eps-regression", epsilon = eps, 
                       cost = input$svr_cost, kernel = "linear", scale = FALSE)
        fitted_vals <- predict(svr_fit, newdata = x)
        calc_pred_obs_slope(y, fitted_vals)
      }, error = function(e) NA)
    })
    
    r_squared <- sapply(eps_values, function(eps) {
      tryCatch({
        svr_fit <- svm(x, y, type = "eps-regression", epsilon = eps, 
                       cost = input$svr_cost, kernel = "linear", scale = FALSE)
        fitted_vals <- predict(svr_fit, newdata = x)
        calc_r_squared(y, fitted_vals)
      }, error = function(e) NA)
    })
    
    # Remove NA values
    valid_idx <- !is.na(slopes) & !is.na(r_squared)
    slope_data <- data.frame(
      epsilon = eps_values[valid_idx], 
      slope = slopes[valid_idx],
      r_squared = r_squared[valid_idx]
    )
    
    if (nrow(slope_data) > 0) {
      # Scale R² to fit nicely with slope
      r2_scaled <- slope_data$r_squared * (max(slope_data$slope) - min(slope_data$slope)) / 
                   (max(slope_data$r_squared) - min(slope_data$r_squared)) + min(slope_data$slope)
      
      ggplot(slope_data, aes(x = epsilon)) +
        geom_line(aes(y = slope, color = "Slope"), linewidth = 1.5) +
        geom_line(aes(y = r2_scaled, color = "R²"), linewidth = 1.5, linetype = "dashed") +
        geom_point(aes(y = slope, color = "Slope"), size = 1.5, alpha = 0.7) +
        geom_point(aes(y = r2_scaled, color = "R²"), size = 1.5, alpha = 0.7) +
        geom_hline(yintercept = 1, linetype = "dotted", color = "red", alpha = 0.7) +
        geom_vline(xintercept = input$svr_epsilon, linetype = "dotted", color = "blue") +
        scale_color_manual(values = c("Slope" = "#6f42c1", "R²" = "#17a2b8")) +
        scale_y_continuous(
          name = "Slope (Pred vs Obs)",
          sec.axis = sec_axis(
            trans = ~ (. - min(slope_data$slope)) / (max(slope_data$slope) - min(slope_data$slope)) * 
                      (max(slope_data$r_squared) - min(slope_data$r_squared)) + min(slope_data$r_squared),
            name = "R² (Goodness of Fit)"
          )
        ) +
        labs(title = paste0("SVR: Slope & R² vs ε (Cost=", input$svr_cost, ")"),
             x = "SVR Epsilon ε", color = "Metric") +
        theme_minimal() +
        theme(legend.position = "bottom")
    } else {
      ggplot() + 
        labs(title = "SVR: No valid results", x = "SVR Epsilon ε", y = "Slope") +
        theme_minimal()
    }
  })
  
  # Data table
  output$data_table <- DT::renderDataTable({
    data <- bmr_data()
    display_data <- data[, .(
      Weight = round(weight, 1),
      BMR_True = round(bmr_true, 1),
      BMR_Observed = round(bmr_obs, 1),
      Error = round(bmr_obs - bmr_true, 1),
      Outlier = ifelse(is_outlier, "Yes", "No")
    )]
    
    # Color outliers in the table
    datatable(display_data, options = list(pageLength = 10, scrollX = TRUE)) %>%
      formatStyle("Outlier",
                  backgroundColor = styleEqual("Yes", "rgba(255, 0, 0, 0.1)"))
  }, options = list(pageLength = 10, scrollX = TRUE))
}

# Run the app
shinyApp(ui = ui, server = server)
