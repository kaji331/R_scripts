library(shiny)
library(bs4Dash)

shinyApp(
    ui=bs4DashPage(
        header=bs4DashNavbar(title="test"),
        sidebar=bs4DashSidebar(
            menuItem(
                "opts",
                sliderInput("x",
                            "x",
                            min=0,
                            max=10,
                            value=0,
                            step=1)
            ),
            skin="light"
        ),
        body=bs4DashBody(
            plotOutput("y")
        )
    ),

server=function(input,output,session) {
                output$y <- renderPlot({
                    plot(0:input$x,col="red")
                })
})
