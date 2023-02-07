currentTime <- function(id) {
    ns <- NS(id)
    textOutput(ns("ct"))
}

currentTimeServer <- function(id) {
    moduleServer(id,
                 function(input,output,session) {
                     output$ct <- reactive({
                         Sys.time()
                     })
                 })
}
