############################################################################################
# Shiny app for simple SEIR model implementation
#
# Miqdad Asaria
# April 2020
############################################################################################

library(shiny)
library(DT)

source("SEIR.R")

ui = fluidPage(
    titlePanel("Simple SEIR model"),

    sidebarLayout(
        sidebarPanel(
            numericInput("N",
                         "Population (N):",
                         value = 100000, 
                         min = 1, 
                         max = NA, 
                         step = 1),
            numericInput("R0",
                         "Basic reproductive number (R0):",
                         min = 0,
                         max = 20,
                         value = 2.2, 
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
            numericInput("delta_t",
                        "Time step in the model (dt):",
                        min = 0,
                        max = 5,
                        value = 0.1,
                        step = 0.1),
            numericInput("days",
                        "Number of days to run model for:",
                        min = 10,
                        max = 1000,
                        value = 150,
                        step = 10),
            tags$div(
                HTML("<small><small>
         <p>This site was produced by <a href='http://www.lse.ac.uk/lse-health/people/miqdad-asaria'>Miqdad Asaria</a> 
         <p>Source code can be found <a href='https://github.com/miqdadasaria/infectious_disease_modelling'>here</a>.
         </small></small>")
            )
        ),

        mainPanel(
            tabsetPanel(id="tabset",
                        tabPanel("SEIR Plot", plotOutput("seir_plot"), 
                                 htmlOutput("summary")),
                        tabPanel("Markov Trace", div(dataTableOutput("seir_markov_trace"), style = "font-size:70%"))
            )
        )
    )
)

# Define server logic required to draw a histogram
server = function(input, output) {
    seirModel = reactive({
        simulate_SEIR(input$N, 
                     input$R0, 
                     input$D, 
                     input$D_pre, 
                     input$delta_t, 
                     input$days)
    })
    
    output$seir_plot = renderPlot({
        seirModel()$seir_plot
    })
    
    output$seir_markov_trace = renderDataTable({
        datatable(seirModel()$seir_markov_trace,
                  style = 'bootstrap',
                  rownames = FALSE,
                  options = list(pageLength = 20, autoWidth = TRUE, dom='ftrpi'))
    })
    
    output$summary = renderText({
        data=tail(seirModel()$seir_markov_trace,1)
        HIT = seirModel()$HIT
        HIT_date = seirModel()$HIT_date
        paste0("After ", round(data$time), " days:<br>",
               "<font color='blue'><b>Infected overall: </b></font> ", round(data$Total_I)," (<b>",round(100*data$Total_I/input$N),"%</b>)<br>", 
               "<font color='green'><b>Recovered overall: </b></font> ", round(data$R), " (<b>",round(100*data$R/input$N),"%</b>)<br>", 
               "<font color='red'><b>Susceptible at end: </b></font> ", round(data$S), " (<b>",round(100*data$S/input$N),"%</b>)<br>",
               "The <b>herd immunity threshold</b> is: <b>", round(100*HIT),"%</b>", 
               if(is.na(HIT_date)){" not achieved in this simulation"} else {paste0(" achieved after <b>", HIT_date, "</b> days after which number of new infections fall")}
        )
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
