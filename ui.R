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
      
      # Input: Interval, number of clusters
      sliderInput("nclusters", "Number of clusters:",
                  min = 6, max = 500,
                  value = 6, step = 1),
      
      # Input: Interval, number of subjects per cluster-period
      sliderInput("nsubjects", "Number of subjects per cluster-period:",
                  min = 5, max = 150,
                  value = 50, step = 5),
      
      # Input: Interval, base correlation rho0
      sliderInput("rho0", "Base correlation:",
                  min = 0, max = 0.2,
                  value = 0.023, step = 0.001),
      
      # Input: Checkbox for whether to view power
      checkboxInput("viewpower", "View power?", FALSE),
      
      conditionalPanel(
        condition = "input.viewpower == true",
        # Input: Interval, effect size for power calculation
        sliderInput("effsize", "Effect size:",
                    min = 0.01, max = 1.0,
                    value = 0.2, step = 0.01)
      ),
      
      # Update button to defer the rendering of output until user
      # clicks the button
      actionButton("update", "Update")
    ),
    
    # Main panel for displaying output
    mainPanel(
      
      # Output: Plot of variance under CT model
      
      uiOutput("plotheader1a"),
      
      uiOutput("plotheader1b"),
      
      plotlyOutput("plot1"),
      
      # Output: Plot of relative variance, HH vs CT vs models

      uiOutput("plotheader2a"),
      
      uiOutput("plotheader2b"),

      plotlyOutput("plot2"),
      
      # Output: Plot of relative variance, DT vs CT models

      uiOutput("plotheader3a"),
      
      uiOutput("plotheader3b"),

      plotlyOutput("plot3"),
      
      # Output: Plot of power under CT model
      
      uiOutput("plotheader4a"),
      
      uiOutput("plotheader4b"),
      
      plotlyOutput("plot4")
    )
  )
))
