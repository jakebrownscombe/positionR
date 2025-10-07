library(shiny)
library(plotly)
library(raster)
library(DT)
library(shinydashboard)

# Source UI and server components
source("ui.R")
source("server.R")

# Run the application
shinyApp(ui = ui, server = server)
