library(shinyFiles)
source('omicsViewer.R')
folder_dir <- "../appData"
server <- function(input, output, session) {
    imputationMethod = 'chen-meng'
    if (file.exists(file.path(folder_dir,'config.yaml'))) {
    config <- yaml::read_yaml(file.path(folder_dir,'config.yaml'))
    imputationMethod <- config$normalization$imputationMethod
    print(paste0('imputation method was chaned to :', imputationMethod ))
    } else {
    print("Config File not found.")
    }
    callModule(app_module, id = "app", .dir = reactive(folder_dir), esetLoader = readESVObj,imputationMethod = imputationMethod)
 
}


