# User interface for continuous-time decay Shiny app
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

library(shiny)
library(plotly)

shinyUI(fluidPage(
  
  # App title
  titlePanel("Continuous-time correlation decay"),
  
  # Sidebar layout
  sidebarLayout(
    
    # Sidebar panel for inputs
    sidebarPanel(
      
      # Input: Interval, number of time periods
      sliderInput("nperiods", "Number of periods:",
                  min = 3, max = 8,
                  value = 4, step = 1,
                  ticks=FALSE),
      
      # Input: Interval, number of subjects per cluster-period
      sliderInput("nsubjects", "Number of subjects per cluster-period:",
                  min = 5, max = 250,
                  value = 50, step = 5),
      
      # Input: Interval, base correlation rho0
      sliderInput("rho0", "Base correlation:",
                  min = 0, max = 0.2,
                  value = 0.023, step = 0.001),
      
      # Update button to defer the rendering of output until user
      # clicks the button
      actionButton("update", "Update")
    ),
    
    # Main panel for displaying output
    mainPanel(
      
      # Output: Plotly plot of variance under CT model
      
      uiOutput("plotheader1a"),
      
      uiOutput("plotheader1b"),
      
      plotlyOutput("plot1"),
      
      # Output: Plotly plot of relative variance, HH vs CT vs models

      uiOutput("plotheader2a"),
      
      uiOutput("plotheader2b"),

      plotlyOutput("plot2"),
      
      # Output: Plotly plot of relative variance, DT vs CT models

      uiOutput("plotheader3a"),
      
      uiOutput("plotheader3b"),

      plotlyOutput("plot3")
    )
  )
))
