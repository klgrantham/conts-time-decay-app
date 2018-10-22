# Underlying calculations and plot generation for
# continuous-time decay Shiny app
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

library(shiny)
library(plotly)

source('conts_time_decay.R')

shinyServer(function(input, output, session) {

  observe({
    updateNumericInput(session, inputId="nclusters", value=2*(input$nperiods - 1))
  })
  
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
    
    res <- generate_var_results_prog(input$nperiods, input$nclusters,
                                     input$nsubjects, input$rho0, updateProgress)
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
    tags$h4("Continuous-time correlation decay (CCD)")
  })
  
  output$plot1 <- renderPlotly({
    res <- getresults()
    p <- plot_ly(res, height=500, width=800, x=~decay, y=~ctpllel, name="Parallel", type="scatter",
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
    tags$h4("Uniform correlation (UC) vs continuous-time correlation decay (CCD)")
  })
  
  output$plot2 <- renderPlotly({
    res <- getresults()
    relres <- res %>% 
      mutate(ratioSW=HHSW/ctSW, ratiocrxo=HHcrxo/ctcrxo, ratiopllel=HHpllel/ctpllel) %>%
      select(decay, rhoUC, starts_with('ratio'))
    p <- plot_ly(relres, height=500, width=800, x=~decay, y=~ratiopllel, name="Parallel", type="scatter",
                 mode="lines", hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=16)),
                 text=~paste("Decay:", round(decay, 2), "<br>Relative variance:", round(ratiopllel, 3),
                             "<br>rhoUC:", round(rhoUC, 3)),
                 line=list(color="#00BA38", width=4, dash="dash")) %>%
      add_trace(y=~ratiocrxo, name="CRXO", hoverinfo="text",
                text=~paste("Decay:",  round(decay, 2), "<br>Relative variance:", round(ratiocrxo, 3),
                            "<br>rhoUC:", round(rhoUC, 3)),
                line=list(color="#F8766D", width=4, dash="dashdot")) %>%
      add_trace(y=~ratioSW, name="SW", hoverinfo="text",
                text=~paste("Decay:", round(decay, 2), "<br>Relative variance:", round(ratioSW, 3),
                            "<br>rhoUC:", round(rhoUC, 3)),
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
    tags$h4("Discrete-time vs continuous-time correlation decay")
  })
  
  output$plot3 <- renderPlotly({
    res <- getresults()
    relres <- res %>% 
      mutate(ratioSW=dtSW/ctSW, ratiocrxo=dtcrxo/ctcrxo, ratiopllel=dtpllel/ctpllel) %>%
      select(decay, starts_with('ratio'))
    p <- plot_ly(relres, height=500, width=800, x=~decay, y=~ratiopllel, name="Parallel", type="scatter",
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
  
  output$plotheader4a <- eventReactive(input$update, {
    header4a()
  })
  
  output$plotheader4b <- eventReactive(input$update, {
    header4b()
  })
  
  header4a <- renderPrint({
    tags$h3(paste0("Power to detect effect size of ", input$effsize))
  })
  
  header4b <- renderPrint({
    tags$h4("Continuous-time correlation decay (CCD)")
  })
  
  output$plot4 <- renderPlotly({
    res <- getresults()
    siglevel <- 0.05
    pow <- powdf(res, input$effsize, siglevel)
    p <- plot_ly(pow, height=500, width=800, x=~decay, y=~ctpllel, name="Parallel", type="scatter",
                 mode="lines", hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=16)),
                 text=~paste("Decay:", round(decay, 2), "<br>Power:", round(ctpllel, 3)),
                 line=list(color="#00BA38", width=4, dash="dash")) %>%
      add_trace(y=~ctcrxo, name="CRXO", hoverinfo="text",
                text=~paste("Decay:",  round(decay, 2), "<br>Power:", round(ctcrxo, 3)),
                line=list(color="#F8766D", width=4, dash="dashdot")) %>%
      add_trace(y=~ctSW, name="SW", hoverinfo="text",
                text=~paste("Decay:", round(decay, 2), "<br>Power:", round(ctSW, 3)),
                line=list(color="#619CFF", width=4, dash="solid")) %>%
      add_trace(y=1, line=list(color="black", width=1, dash="solid"), hoverinfo="none", showlegend=FALSE) %>%
      layout(xaxis=list(title="Decay", titlefont=list(size=18), showline=TRUE,
                        tickmode="auto", tickfont=list(size=16), nticks=6, ticks="inside",
                        mirror=TRUE, showgrid=FALSE),
             yaxis=list(title="Power", titlefont=list(size=18), tickfont=list(size=16),
                        mirror=TRUE, showline=TRUE),
             legend=list(orientation="h", xanchor="center", x=0.5, y=-0.3, font=list(size=16)),
             margin=list(l=100, r=40))
    p$elementId <- NULL
    p
  })
  
  getresultstable <- eventReactive(input$update, {
    res <- getresults()
    siglevel <- 0.05
    pow <- powdf(res, input$effsize, siglevel)
    results <- data.frame(decay=res$decay, rhoCCD=signif(input$rho0, 3), rhoUC=signif(res$rhoUC, 3),
                          CCD_crxo=signif(res$ctcrxo, 3), CCD_pllel=signif(res$ctpllel, 3),
                          CCD_SW=signif(res$ctSW, 3), UC_crxo=signif(res$HHcrxo, 3),
                          UC_pllel=signif(res$HHpllel, 3), UC_SW=signif(res$HHSW, 3),
                          power_CCD_crxo=signif(pow$ctcrxo, 3), power_CCD_pllel=signif(pow$ctpllel, 3),
                          power_CCD_SW=signif(pow$ctSW, 3)
                )
  })

  output$resultstable <- renderTable({
    restable <- getresultstable()
    restable$rhoCCD <- format(restable$rhoCCD, 3)
    restable$rhoUC <- format(restable$rhoUC, 3)
    restable$CCD_crxo <- format(restable$CCD_crxo, 3)
    restable$CCD_pllel <- format(restable$CCD_pllel, 3)
    restable$CCD_SW <- format(restable$CCD_SW, 3)
    restable$UC_crxo <- format(restable$UC_crxo, 3)
    restable$UC_pllel <- format(restable$UC_pllel, 3)
    restable$UC_SW <- format(restable$UC_SW, 3)
    restable$power_CCD_crxo <- format(restable$power_CCD_crxo, 3)
    restable$power_CCD_pllel <- format(restable$power_CCD_pllel, 3)
    restable$power_CCD_SW <- format(restable$power_CCD_SW, 3)
    names(restable) <- c("Decay", "rhoCCD", "rhoUC", "CCD CRXO", "CCD Parallel", "CCD SW",
                         "UC CRXO", "UC Parallel", "UC SW", "Power (CCD CRXO)",
                         "Power (CCD Parallel)", "Power (CCD SW)")
    restable[order(restable$Decay),]
  })
  
  output$resultsdownload <- downloadHandler(
    filename = function(){
      if(input$rho0==1){rho0char <- 100}else{rho0char <- strsplit(as.character(input$rho0),"\\.")[[1]][2]}
      if(input$effsize==1){effsizechar <- 100}else{effsizechar <- strsplit(as.character(input$effsize),"\\.")[[1]][2]} 
      paste0("results_", "N", input$nperiods, "T", input$nperiods, "m",
             input$nsubjects, "rho", rho0char, "effsize", effsizechar, ".csv")
    },
    content = function(file){
      write.csv(getresultstable(), file, row.names=TRUE)
    }
  )

  # Tables showing design matrices
  
  # Retrieve design matrices
  getdesmats <- eventReactive(input$update, {
    crxo <- crxodesmat(input$nperiods, input$nclusters)
    pllel <- plleldesmat(input$nperiods, input$nclusters)
    SW <- SWdesmat(input$nperiods, input$nclusters)
    list(crxo, pllel, SW)
  })
  
  # Retrieve design matrices
  getdesmats <- function(){
    crxo <- crxodesmat(input$nperiods, input$nclusters)
    pllel <- plleldesmat(input$nperiods, input$nclusters)
    SW <- SWdesmat(input$nperiods, input$nclusters)
    list(crxo, pllel, SW)
  }
  

  output$crxodesmat <- renderTable({
    desmats <- getdesmats()
    desmats[[1]]
    },
    colnames=FALSE, digits=0, spacing='xs'
  )

  output$plleldesmat <- renderTable({
    desmats <- getdesmats()
    desmats[[2]]
    },
    colnames=FALSE, digits=0, spacing='xs'
  )
  
  output$SWdesmat <- renderTable({
    desmats <- getdesmats()
    desmats[[3]]
    },
    colnames=FALSE, digits=0, spacing='xs'
  )

  output$plotheader5a <- renderPrint({
    tags$h4("Cluster randomised crossover (CRXO)")
  })
  
  output$plotheader5b <- renderPrint({
    tags$h4("Parallel")
  })
  
  output$plotheader5c <- renderPrint({
    tags$h4("Stepped wedge (SW)")
  })
})
