# Underlying calculations and plot generation for
# continuous-time decay Shiny app
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

library(shiny)
library(plotly)

source('conts_time_decay.R')

shinyServer(function(input, output) {
  
  getresults <- eventReactive(input$update, {
    
    # Create a Progress object
    progress <- shiny::Progress$new(style = "old")
    progress$set(message = "Generating results...", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    # Create a closure to update progress.
    # Each time this is called:
    # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
    #   distance. If non-NULL, it will set the progress to that value.
    # - It also accepts optional detail text.
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 5
      }
      progress$set(value = value, detail = detail)
    }
    
    res <- generate_var_results_prog(input$nperiods, input$nsubjects, input$rho0, updateProgress)
    return(res)
  })
  
  output$plotheader1a <- eventReactive(input$update, {
    header1a()
  })
  
  output$plotheader1b <- eventReactive(input$update, {
    header1b()
  })
  
  header1a <- renderPrint({
    tags$h3("Variance of treatment effect estimator")
  })
  
  header1b <- renderPrint({
    tags$h4("Continuous-time decaying correlation")
  })
  
  output$plot1 <- renderPlotly({
    res <- getresults()
    p <- plot_ly(res, height=400, width=650, x=~decay, y=~ctpllel, name="Parallel", type="scatter",
                 mode="lines", hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=16)),
                 text=~paste("Decay:", round(decay, 2), "<br>Variance:", round(ctpllel, 3)),
                 line=list(color="#00BA38", width=4, dash="dash")) %>%
      add_trace(y=~ctcrxo, name="CRXO", hoverinfo="text",
                text=~paste("Decay:", round(decay, 2), "<br>Variance:", round(ctcrxo, 3)),
                line=list(color="#F8766D", width=4, dash="dashdot")) %>%
      add_trace(y=~ctSW, name="SW", hoverinfo="text",
                text=~paste("Decay:", round(decay, 2), "<br>Variance:", round(ctSW, 3)),
                line=list(color="#619CFF", width=4, dash="solid")) %>%
      layout(xaxis=list(title="Decay", titlefont=list(size=18), showline=TRUE,
                        tickmode="auto", tickfont=list(size=16), nticks=6, ticks="inside",
                        mirror=TRUE, showgrid=FALSE),
             yaxis=list(title="Variance", titlefont=list(size=18), tickfont=list(size=16),
                        mirror=TRUE, showline=TRUE),
             legend=list(orientation="h", xanchor="center", yanchor="bottom", x=0.5, y=-0.5, font=list(size=16)),
             margin=list(l=100, r=40))
    p$elementId <- NULL # Workaround to suppress warning due to an incompatility between shiny and plotly
    p
  })
  
  output$plotheader2a <- eventReactive(input$update, {
    header2a()
  })
  
  output$plotheader2b <- eventReactive(input$update, {
    header2b()
  })
  
  header2a <- renderPrint({
    tags$h3("Relative variance of treatment effect estimator")
  })
  
  header2b <- renderPrint({
    tags$h4("Uniform vs Continuous-time decaying correlation")
  })
  
  output$plot2 <- renderPlotly({
    res <- getresults()
    relres <- res %>% 
      mutate(ratioSW=HHSW/ctSW, ratiocrxo=HHcrxo/ctcrxo, ratiopllel=HHpllel/ctpllel) %>%
      select(decay, starts_with('ratio'))
    p <- plot_ly(relres, height=400, width=650, x=~decay, y=~ratiopllel, name="Parallel", type="scatter",
                 mode="lines", hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=16)),
                 text=~paste("Decay:", round(decay, 2), "<br>Relative variance:", round(ratiopllel, 3)),
                 line=list(color="#00BA38", width=4, dash="dash")) %>%
      add_trace(y=~ratiocrxo, name="CRXO", hoverinfo="text",
                text=~paste("Decay:",  round(decay, 2), "<br>Relative variance:", round(ratiocrxo, 3)),
                line=list(color="#F8766D", width=4, dash="dashdot")) %>%
      add_trace(y=~ratioSW, name="SW", hoverinfo="text",
                text=~paste("Decay:", round(decay, 2), "<br>Relative variance:", round(ratioSW, 3)),
                line=list(color="#619CFF", width=4, dash="solid")) %>%
      add_trace(y=1, line=list(color="black", width=1, dash="solid"), hoverinfo="none", showlegend=FALSE) %>%
      layout(xaxis=list(title="Decay", titlefont=list(size=18), showline=TRUE,
                        tickmode="auto", tickfont=list(size=16), nticks=6, ticks="inside",
                        mirror=TRUE, showgrid=FALSE),
             yaxis=list(title="Relative variance", titlefont=list(size=18), tickfont=list(size=16),
                        mirror=TRUE, showline=TRUE),
             legend=list(orientation="h", xanchor="center", x=0.5, y=-0.3, font=list(size=16)),
             margin=list(l=100, r=40))
    p$elementId <- NULL
    p
  })
  
  output$plotheader3a <- eventReactive(input$update, {
    header3a()
  })
  
  output$plotheader3b <- eventReactive(input$update, {
    header3b()
  })
  
  header3a <- renderPrint({
    tags$h3("Relative variance of treatment effect estimator")
  })
  
  header3b <- renderPrint({
    tags$h4("Discrete-time vs Continuous-time decaying correlation")
  })
  
  output$plot3 <- renderPlotly({
    res <- getresults()
    relres <- res %>% 
      mutate(ratioSW=dtSW/ctSW, ratiocrxo=dtcrxo/ctcrxo, ratiopllel=dtpllel/ctpllel) %>%
      select(decay, starts_with('ratio'))
    p <- plot_ly(relres, height=400, width=650, x=~decay, y=~ratiopllel, name="Parallel", type="scatter",
                 mode="lines", hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=16)),
                 text=~paste("Decay:", round(decay, 2), "<br>Relative variance:", round(ratiopllel, 3)),
                 line=list(color="#00BA38", width=4, dash="dash")) %>%
      add_trace(y=~ratiocrxo, name="CRXO", hoverinfo="text",
                text=~paste("Decay:",  round(decay, 2), "<br>Relative variance:", round(ratiocrxo, 3)),
                line=list(color="#F8766D", width=4, dash="dashdot")) %>%
      add_trace(y=~ratioSW, name="SW", hoverinfo="text",
                text=~paste("Decay:", round(decay, 2), "<br>Relative variance:", round(ratioSW, 3)),
                line=list(color="#619CFF", width=4, dash="solid")) %>%
      add_trace(y=1, line=list(color="black", width=1, dash="solid"), hoverinfo="none", showlegend=FALSE) %>%
      layout(xaxis=list(title="Decay", titlefont=list(size=18), showline=TRUE,
                        tickmode="auto", tickfont=list(size=16), nticks=6, ticks="inside",
                        mirror=TRUE, showgrid=FALSE),
             yaxis=list(title="Relative variance", titlefont=list(size=18), tickfont=list(size=16),
                        mirror=TRUE, showline=TRUE),
             legend=list(orientation="h", xanchor="center", x=0.5, y=-0.3, font=list(size=16)),
             margin=list(l=100, r=40))
    p$elementId <- NULL
    p
  })
})
