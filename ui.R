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
      sliderInput("nperiods", "Number of periods, T:",
                  min = 3, max = 8,
                  value = 4, step = 1,
                  ticks=FALSE),
      
      # Input: Interval, number of clusters
      numericInput("nclusters", "Number of clusters, N:",
                  min = 6, max = 500,
                  value = 6, step = 2),
      
      # Input: Interval, number of subjects per cluster-period
      sliderInput("nsubjects", "Number of subjects per cluster-period, m:",
                  min = 5, max = 150,
                  value = 50, step = 5),
      
      # Input: Interval, base correlation rho0
      sliderInput("rho0",
                  label = HTML(paste0("Base correlation, &rho;", ":")),
                  min = 0, max = 0.2,
                  value = 0.023, step = 0.001),
      
      # Input: Interval, effect size for power calculation
      sliderInput("effsize", "Effect size:",
                  min = 0.05, max = 1.0,
                  value = 0.2, step = 0.05),
      
      # Update button to defer the rendering of output until user
      # clicks the button
      actionButton("update", "Update")
    ),
    
    # Main panel for displaying output
    mainPanel(
      
      tabsetPanel(type = "tabs",
                  
                  tabPanel("Variance",
                           uiOutput("plotheader1a"), uiOutput("plotheader1b"),
                           plotlyOutput("plot1")),
                  
                  tabPanel("Relative variance",
                           uiOutput("plotheader2a"), uiOutput("plotheader2b"),
                           plotlyOutput("plot2")),
                  
                  tabPanel("Power",
                           uiOutput("plotheader4a"), uiOutput("plotheader4b"),
                           plotlyOutput("plot4")),
                  
                  tabPanel("Results table",
                           br(),
                           downloadButton("resultsdownload", "Download results (.csv)"),
                           hr(),
                           tableOutput("resultstable")),
                  
                  tabPanel("Design matrices",
                           uiOutput("plotheader5a"), tableOutput("crxodesmat"),
                           uiOutput("plotheader5b"), tableOutput("plleldesmat"),
                           uiOutput("plotheader5c"), tableOutput("SWdesmat"))
      )
    )
  )
))
