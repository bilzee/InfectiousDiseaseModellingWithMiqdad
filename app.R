############################################################################################
# Shiny app for simple SEIR model implementation
#
# Miqdad Asaria
# April 2020
############################################################################################

library("shiny")
library("DT")
library("plotly")

source("SEIR.R")
source("test.R")
source("formatting_functions.R")

ui = fluidPage(
    titlePanel("Simple SEIR model"),

    sidebarLayout(
        sidebarPanel(
            numericInput("N",
                         "Population (N):",
                         value = 60000000, 
                         min = 1, 
                         max = NA, 
                         step = 1),
            numericInput("R0",
                         HTML("Basic reproductive number (R<sub>0</sub>):"),
                         min = 0,
                         max = 20,
                         value = 3.25, 
                         step = 0.1),
            numericInput("D_pre",
                        "Average duratiuon of pre-infectious period (D'):",
                        min = 0,
                        max = 100,
                        value = 5.2,
                        step = 0.1),
            numericInput("D",
                        "Average duratiuon of infectious period (D):",
                        min = 0,
                        max = 100,
                        value = 2.9,
                        step = 0.1),
            numericInput("R0_sup",
                         HTML("Policy adjusted reproductive number (R<sub>sup</sub>):"),
                         min = 0,
                         max = 20,
                         value = 0.6, 
                         step = 0.1),
            numericInput("sup_start",
                         HTML("Policy start (days from start of epidemic):"),
                         min = 0,
                         max = 1000,
                         value = 80, 
                         step = 1),
            numericInput("D_sup",
                         HTML("Policy duration (days):"),
                         min = 0,
                         max = 1000,
                         value = 60, 
                         step = 1),
            numericInput("days",
                         "Number of days to run model for:",
                         min = 10,
                         max = 1000,
                         value = 365,
                         step = 10),
            tags$hr(),
            sliderInput("sensitivity",
                        "Sensitivity of test (%)",
                        min = 0,
                        max = 100,
                        value = 93.8,
                        step = 0.1),
            sliderInput("specificity",
                        "Specificity of test (%)",
                        min = 0,
                        max = 100,
                        value = 95.6,
                        step = 0.1),
            sliderInput("max_prevalence",
                        "Max prevalence for test plot (%)",
                        min = 0,
                        max = 100,
                        value = 10.0,
                        step = 0.1),
            #tags$hr,
            tags$div(
                HTML("<small><small>
         <p>This site was produced by <a href='http://www.lse.ac.uk/lse-health/people/miqdad-asaria'>Miqdad Asaria</a> 
         <p>Source code can be found <a href='https://github.com/miqdadasaria/infectious_disease_modelling'>here</a>.
         </small></small>")
            )
        ),

        mainPanel(
            tabsetPanel(id="tabset",
                        tabPanel("The \"Curve\"", plotlyOutput("infection_plot")),
                        tabPanel("SEIR Plot (supression)", plotlyOutput("seir_sup_plot"), 
                                 htmlOutput("sup_summary")),                        
                        #tabPanel("Markov Trace (supression)", div(dataTableOutput("seir_sup_markov_trace"), style = "font-size:70%")),
                        tabPanel("SEIR Plot (do nothing)", plotlyOutput("seir_plot"), 
                                 htmlOutput("summary")),
                        #tabPanel("Markov Trace (do nothing)", div(dataTableOutput("seir_markov_trace"), style = "font-size:70%")),
                        tabPanel("Tests", plotlyOutput("test_plot")),
                        tabPanel("HIT", plotlyOutput("hit_plot"))
            )
        )
    )
)

# Define server logic required to draw a histogram
server = function(input, output) {
    seir_model = reactive({
        simulate_SEIR_differential(input$N, 
                     input$R0, 
                     input$D, 
                     input$D_pre, 
                     input$days)
    })
    
    seir_sup_model = reactive({
        supress_release_SEIR(N = input$N, 
                             R0 = input$R0, 
                             R0_sup = input$R0_sup, 
                             D = input$D, 
                             D_pre = input$D_pre, 
                             D_sup = input$D_sup, 
                             sup_start = input$sup_start, 
                             total_days = input$days)
    })
    
    output$seir_plot = renderPlotly({
        convert_ggplot_to_plotly(plot_SEIR(seir_model(),
                           input$N, 
                           input$R0, 
                           input$D, 
                           input$D_pre))
    })
    
    output$seir_sup_plot = renderPlotly({
        convert_ggplot_to_plotly(plot_SEIR(seir_sup_model(),
                           input$N, 
                           input$R0, 
                           input$D, 
                           input$D_pre))
    })
    
    output$infection_plot = renderPlotly({
        convert_ggplot_to_plotly(plot_infection_curves(seir_model(), 
                                             seir_sup_model(),
                                             paste0("Baseline R0 = ",input$R0),
                                             paste0("Supression R0 reduced to ", input$R0_sup, " from day ", input$sup_start, " for ", input$D_sup, " days")))
    })
    
    output$seir_markov_trace = renderDataTable({
        markov_trace = round(seir_model(),2)
        datatable(markov_trace,
                  style = 'bootstrap',
                  rownames = FALSE,
                  options = list(pageLength = 20, autoWidth = TRUE, dom='ftrpi'))
    })
    
    output$seir_sup_markov_trace = renderDataTable({
        markov_trace = round(seir_sup_model(),2)
        datatable(markov_trace,
                  style = 'bootstrap',
                  rownames = FALSE,
                  options = list(pageLength = 20, autoWidth = TRUE, dom='ftrpi'))
    })
    
    output$summary = renderText({
        data=tail(seir_model(),1)
        HIT = calculate_HIT(input$R0)
        HIT_date = calculate_HIT_date(seir_model(),input$R0,input$N)
        f = 1 / input$D_pre
        r = 1 / input$D
        beta = (input$R0 * r) / input$N 
        paste0("Key derived model parameters: <b>\u03B2</b> = <b>",round(beta,10),"</b>, <b>f</b> = <b>",round(f,5),"</b>, <b>r</b> = <b>", round(r,5), "</b><br><br>",
            "After ", round(data$time), " days:<br>",
               "<font color='blue'><b>Infected overall: </b></font> ", round(data$Total_I)," (<b>",round(100*data$Total_I/input$N),"%</b>)<br>", 
               "<font color='green'><b>Recovered overall: </b></font> ", round(data$R), " (<b>",round(100*data$R/input$N),"%</b>)<br>", 
               "<font color='red'><b>Susceptible at end: </b></font> ", round(data$S), " (<b>",round(100*data$S/input$N),"%</b>)<br>",
               "The <b>herd immunity threshold</b> is: <b>", round(100*HIT),"%</b>", 
               if(is.na(HIT_date)){" not achieved in this simulation"} else {paste0(" achieved after <b>", HIT_date, "</b> days after which number of new infections fall")}
        )
    })
    
    
    output$sup_summary = renderText({
        data=tail(seir_sup_model(),1)
        HIT = calculate_HIT(input$R0)
        HIT_date = calculate_HIT_date(seir_sup_model(),input$R0,input$N)
        paste0("After ", round(data$time), " days:<br>",
               "<font color='blue'><b>Infected overall: </b></font> ", round(data$Total_I)," (<b>",round(100*data$Total_I/input$N),"%</b>)<br>", 
               "<font color='green'><b>Recovered overall: </b></font> ", round(data$R), " (<b>",round(100*data$R/input$N),"%</b>)<br>", 
               "<font color='red'><b>Susceptible at end: </b></font> ", round(data$S), " (<b>",round(100*data$S/input$N),"%</b>)<br>",
               "The <b>herd immunity threshold</b> at unsupressed R<sub>0</sub> of ",input$R0," is: <b>", round(100*HIT),"%</b>", 
               if(is.na(HIT_date)){" not achieved in this simulation"} else {paste0(" achieved after <b>", HIT_date, "</b> days after which number of new infections fall")}
        )
    })
    
    output$hit_plot = renderPlotly({
        convert_ggplot_to_plotly(plot_HIT())
    })
    
    output$test_plot = renderPlotly({
        convert_ggplot_to_plotly(calculate_test_ppv_npv(sensitivity = input$sensitivity/100, specificity = input$specificity/100, max_prevalence = input$max_prevalence/100))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
