# UI
ui <- dashboardPage(
  dashboardHeader(title = "BMR Parameter Covariance Explorer", titleWidth = 350),
  
  dashboardSidebar(
    width = 180,
    sidebarMenu(
      menuItem("Analysis", tabName = "analysis", icon = icon("chart-line"))
    ),
    
    # Model Configuration
    h5("Data Generation:"),
    radioButtons("data_method", "Method:",
                 choices = list("Independent" = "independent",
                                "BMI-driven" = "bmi_driven", 
                                "User-specified" = "user_specified"),
                 selected = "independent"),
    
    h5("Sample Options:"),
    numericInput("sample_size", "Sample Size:", value = 100, min = 50, max = 300, step = 25),
    numericInput("noise_level", "Noise %:", value = 5, min = 0, max = 30, step = 1),
    
    radioButtons("gender_subset", "Gender:",
                 choices = list("Both" = "both", "Male" = "male", "Female" = "female"),
                 selected = "female")
  )
  
  , dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side { background-color: #fff; }
        .box { margin-bottom: 5px; }
        .box-body { padding: 10px; }
        .small-box .inner h3 { font-size: 18px; }
        .small-box .inner p { font-size: 12px; }
        .form-group { margin-bottom: 8px; }
      "))
    )
    
    , fluidRow(
      # Data Generation Method Display
      box(title = "Data Generation", status = "primary", solidHeader = TRUE,
          width = 6, height = "580px"
          
          , verbatimTextOutput("data_method_display")
          
          #CL: put these two tables next to each other
          , div(#style = "height: 200px; width: 300px;", 
                h6("Predictor Ranges (Editable)"),
                rHandsontableOutput("ranges_table")),
          div(style = "height: 200px;", 
              h6("Predictor Correlations"),
              rHandsontableOutput("predictor_correlation_table"))
      )
      
      , box(title = "Parameters & Covariances", status = "warning", solidHeader = TRUE,
            width = 6, height = "580px"
            , verbatimTextOutput("current_model", placeholder = TRUE)
            , checkboxGroupInput("centering_vars", "Centering:",
                                 choices = list("Weight" = "weight",
                                                "Height" = "height", 
                                                "Age" = "age"),
                                 selected = NULL, inline = TRUE)
            
            # , div(style = "height: 220px;"
                  , rHandsontableOutput("combined_params_table")#)
            , div(style = "height: 270px; display: flex;",
                  div(style = "width: 50%; padding-right: 10px;",
                      h6("Parameter Correlations"),
                      rHandsontableOutput("param_correlation_table"
                                          # , height = "340px", width = "500px"
                                          )),
                  div(style = "width: 50%; padding-left: 10px;",
                      h6("Parameter Covariances"), 
                      rHandsontableOutput("param_covariance_table"))
            )
      )
      
    )
    
    , fluidRow(
      box(title = "Model Fit", status = "success", solidHeader = TRUE,
          width = 12, height = "200px"
          , fluidRow(
            column(6 
                   , tableOutput("fit_stats")
            ), column(6
                      , verbatimTextOutput("intercept_analysis")
            ))
      )
    )
    
    
    , fluidRow(
      # Residual Plot 
      box(title = "Residuals vs Fitted", status = "success", solidHeader = TRUE,
          width = 3, height = "300px",
          plotOutput("residual_plot", height = "280px")
      ),
      
      # BMR vs Weight  
      box(title = "BMR vs Weight", status = "info", solidHeader = TRUE,
          width = 3, height = "300px",
          plotOutput("plot_weight", height = "230px")
      ),
      
      # BMR vs Height  
      box(title = "BMR vs Height", status = "info", solidHeader = TRUE,
          width = 3, height = "300px",
          plotOutput("plot_height", height = "230px")
      ),
      
      # BMR vs Age
      box(title = "BMR vs Age", status = "info", solidHeader = TRUE,
          width = 3, height = "300px", 
          plotOutput("plot_age", height = "230px")
      )
    )
    
  )
)
