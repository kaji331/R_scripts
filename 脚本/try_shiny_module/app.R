library(shiny)

ui <- fluidPage(currentTime("test_id"))
server <- function(input,output,session) {
    currentTimeServer("test_id")
}

shinyApp(ui,server)
