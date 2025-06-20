.libPaths( file.path(getwd(), "App", "R-portable", "library") )
library(omicsViewer)
shiny::runApp('shiny', launch.browser = TRUE)
