# User interface for continuous-time decay Shiny app
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

library(shiny)
library(plotly)

shinyUI(fluidPage(
  
  # App title
  titlePanel("Continuous-time decaying correlation model"),
  
  # Sidebar layout
  sidebarLayout(
    
    # Sidebar panel for inputs
    sidebarPanel(
      
      # Input: Interval, number of time periods
      sliderInput("nperiods", "Number of periods:",
                  min = 2, max = 10,
                  value = 4, step = 1,
                  ticks=FALSE),
      
      # Input: Interval, number of subjects per cluster-period
      sliderInput("nsubjects", "Number of subjects per cluster-period:",
                  min = 10, max = 200, # Raise max later once working
                  value = 50, step = 10),
      
      # Input: Interval, base correlation rho0
      sliderInput("rho0", "Base correlation:",
                  min = 0, max = 0.2,
                  value = 0.04, step = 0.005),
      
      # Update button to defer the rendering of output until user
      # clicks the button
      actionButton("update", "Update")
    ),
    
    # Main panel for displaying output
    mainPanel(
      
      # Output: Table summarizing the values entered
      tableOutput("values"),
      
      # Output: Plotly plot of variance under CT model
      plotlyOutput("plot1"),
      
      # Output: Plotly plot of relative variance, CT vs HH models
      plotlyOutput("plot2"),
      
      # Output: Plotly plot of relative variance, CT vs DT models
      plotlyOutput("plot3")
    )
  )
))
