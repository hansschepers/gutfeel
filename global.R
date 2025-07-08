library(shiny)
library(shinydashboard)
library(data.table)
library(ggplot2)
library(corrplot)
library(DT)
library(rhandsontable)

# Original Harris-Benedict coefficients (1919)
HARRIS_BENEDICT <- list(
  male = c(intercept = 66.473, weight = 13.7516, height = 5.0033, age = -6.755),
  female = c(intercept = 655.0955, weight = 9.5634, height = 1.8500, age = -4.676)
)

source("functions.R")
