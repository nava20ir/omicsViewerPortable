source('libs.R')


aaFreq <- function (x) 
{
  aa <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
          "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  x <- x[which(nchar(x) > 0)]
  if (length(unique(nchar(x))) > 1) 
    stop("Different length in x!")
  s0 <- strsplit(x, "|")
  s0 <- do.call(rbind, s0)
  fg <- apply(s0, 2, function(x) table(x)[aa])
  fg[is.na(fg)] <- 0
  s <- sweep(fg, 2, colSums(fg, na.rm = TRUE), "/")
  rownames(s) <- aa
  s
}


addHeatmapAnnotation <- function (x, column = TRUE, var.name = "") 
{
  object <- list()
  object$column <- column
  object$var.name <- var.name
  if (is.numeric(x) && is.vector(x) && !is.matrix(x)) {
    object$type <- 1
    object$x <- x
  }
  else if ((is.character(x) && is.vector(x)) || is.factor(x) || 
           is.logical(x)) {
    ll <- value2color(x)
    object$key <- ll
    object$type <- 2
  }
  else if (is.matrix(x) || is.data.frame(x)) {
    ll <- lapply(seq_len(ncol(x)), function(i) value2color(x[, 
                                                             i]))
    cmat <- t(vapply(ll, function(x) x[["color"]], character(length(ll[[1]][["color"]]))))
    keylist <- lapply(ll, "[[", "key")
    mm <- matrix(1, nrow(cmat), ncol(cmat))
    mm <- colCumsums(mm)
    rownames(mm) <- rownames(cmat) <- names(keylist) <- colnames(x)
    object$key <- ll
    object$mm <- mm
    object$cmat <- cmat
    object$type <- 3
  }
  object
}



addHeatmapAnnotation_plot <- function (object, ...) 
{
  x <- object$x
  type <- object$type
  column <- object$column
  var.name <- object$var.name
  ll <- object$key
  mm <- object$mm
  cmat <- object$cmat
  minColShow <- 50
  minRowShow <- 100
  borderColor <- NA
  if (object$type == 2) {
    if (length(object$key$color) < ifelse(column, minColShow, 
                                          minRowShow)) 
      borderColor <- "gray25"
  }
  else if (object$type == 3) {
    if (ncol(cmat) < ifelse(column, minColShow, minRowShow)) 
      borderColor <- "gray25"
  }
  lim <- list(...)
  if (!is.null(lim$xlim)) {
    if (lim$xlim[2] - lim$xlim[1] < minColShow) 
      borderColor <- "gray25"
  }
  if (!is.null(lim$ylim)) {
    if (lim$ylim[2] - lim$ylim[1] < minRowShow) 
      borderColor <- "gray25"
  }
  if (type == 1) {
    if (!column) 
      x <- -x
    barplot(x, horiz = !column, xaxs = "i", yaxs = "i", space = 0, 
            xpd = FALSE, axes = FALSE, ...)
    if (!column) {
      pp <- par("xaxp")
      at <- seq(pp[1], pp[2], length.out = pp[3] + 1)
      lb <- round(rev(seq(-pp[2], -pp[1], length.out = pp[3] + 
                            1)), digits = 3)
      axis(side = 1, at = at, labels = lb, las = 2)
    }
    else {
      axis(side = 4, las = 2)
    }
  }
  else if (type == 2) {
    barplot(rep(1, length(ll$color)), col = ll$color, space = 0, 
            border = borderColor, xaxs = "i", yaxs = "i", axes = FALSE, 
            horiz = !column, xpd = FALSE, ...)
    mtext(side = ifelse(column, 4, 1), at = 0.5, text = var.name, 
          las = 2, line = 0.5)
  }
  else if (type == 3) {
    for (i in nrow(mm):1) barplot(mm[i, ], space = 0, col = cmat[i, 
    ], xaxs = "i", yaxs = "i", axes = FALSE, border = borderColor, 
    add = i != nrow(mm), horiz = !column, xpd = FALSE, 
    ...)
    mtext(side = ifelse(column, 4, 1), at = mm[, 1] - 0.5, 
          text = rownames(mm), las = 2, line = 0.5)
  }
  else stop("Unreal type!")
}


adist <- function (x, method = "pearson") 
{
  arg.dist <- c("euclidean", "maximum", "manhattan", "canberra", 
                "binary", "minkowski")
  arg.cor <- c("spearman", "pearson")
  method <- match.arg(method, choices = c(arg.dist, arg.cor))
  if (method %in% arg.dist) 
    dd <- as.matrix(dist(x, method = method))
  else dd <- 1 - cor(t(x), use = "pairwise")
  med <- median(dd, na.rm = TRUE)
  if (is.na(med)) 
    med <- 0
  dd[is.infinite(dd)] <- med
  dd[is.na(dd)] <- med
  as.dist(dd)
}




app_module <- function (input, output, session, .dir, filePattern = ".(RDS|db|sqlite|sqlite3)$", 
                        additionalTabs = NULL, ESVObj = reactive(NULL), esetLoader = readESVObj, 
                        exprsGetter = getExprs, pDataGetter = getPData, fDataGetter = getFData, 
                        imputeGetter = getExprsImpute, defaultAxisGetter = getAx, 
                        appName = "omicsViewer", appVersion = "1.7.0", imputationMethod = 'custom') 
{
 
  ns <- session$ns
  observe({
    req(.dir())


    
    print(.dir())    
    ll <- list.files(.dir(), pattern = filePattern, ignore.case = TRUE)
    updateSelectizeInput(session = session, inputId = "selectFile", 
                         choices = ll, selected = "")
  })
  reactive_eset <- reactive({
    if (!is.null(ESVObj())) {
      updateSelectizeInput(session, "selectFile", choices = "ESVObj.RDS", 
                           selected = "ESVObj.RDS")
      return(tallGS(ESVObj()))
    }
    req(input$selectFile)
    flink <- file.path(.dir(), input$selectFile)
    sss <- file.size(flink)
    if (sss > 1e+07) 
      show_modal_spinner(text = "Loading data ...")
    v <- esetLoader(flink)
    if (sss > 1e+07) 
      remove_modal_spinner()
    v
  })
  expr <- reactive({
    req(reactive_eset())
    getExprs(reactive_eset())
  })
  pdata <- reactive({
    req(reactive_eset())
    pDataGetter(reactive_eset())
  })
  fdata <- reactive({
    req(reactive_eset())
    fDataGetter(reactive_eset())
  })
  validEset <- function(expr, pd, fd) {
    i1 <- identical(rownames(expr), rownames(fd))
    i2 <- identical(colnames(expr), rownames(pd))
    if (!(i1 && i2)) 
      return(list(FALSE, "The rownames/colnames of exprs not matched to row names of feature data/phenotype data!"))
    TRUE
  }
  print('middle of app module')
  vEset <- reactiveVal(FALSE)
  observe({
    req(expr())
    req(pdata())
    req(fdata())
    x <- validEset(expr = expr(), pd = pdata(), fd = fdata())
    if (!x[[1]]) {
      showModal(modalDialog(title = "Problem in data!", 
                            x[[2]]))
    }
    else {
      vEset(TRUE)
      shinyjs::show("contents")
    }
  })
  d_s_x <- reactive({
    req(eset <- reactive_eset())
    defaultAxisGetter(eset, "sx")
  })
  d_s_y <- reactive({
    req(eset <- reactive_eset())
    defaultAxisGetter(eset, "sy")
  })
  d_f_x <- reactive({
    req(eset <- reactive_eset())
    defaultAxisGetter(eset, "fx")
  })
  d_f_y <- reactive({
    req(eset <- reactive_eset())
    defaultAxisGetter(eset, "fy")
  })
  cormat <- reactive({
    req(eset <- reactive_eset())
    attr(eset, "cormat")
  })
  output$download <- downloadHandler(filename = function() {
    paste0("ExpressenSet", Sys.time(), ".xlsx")
  }, content = function(file) {
    td <- function(tab) {
      ic <- which(vapply(tab, is.list, logical(1)))
      if (length(ic) > 0) {
        for (ii in ic) {
          tab[, ii] <- vapply(tab[, ii], paste, collapse = ";", 
                              FUN.VALUE = character(1))
        }
      }
      id <- rownames(tab)
      if (is.null(id)) 
        id <- paste0("ID", seq_len(nrow(tab)))
      data.frame(ID = id, tab)
    }
    
    ig <- imputeGetter(reactive_eset())
    
    withProgress(message = "Writing table", value = 0, {
      wb <- createWorkbook(creator = "BayBioMS")
      addWorksheet(wb, sheetName = "Phenotype info")
      addWorksheet(wb, sheetName = "Feature info")
      addWorksheet(wb, sheetName = "Expression")
      addWorksheet(wb, sheetName = "Geneset annot")
      incProgress(1/5, detail = "expression matrix")
      writeData(wb, sheet = "Expression", td(expr()))
      if (!is.null(ig)) {
        addWorksheet(wb, sheetName = "Expression_imputed")
        writeData(wb, sheet = "Expression_imputed", td(ig))
      }
      incProgress(1/5, detail = "feature table")
      writeData(wb, sheet = "Feature info", td(fdata()))
      incProgress(1/5, detail = "phenotype table")
      writeData(wb, sheet = "Phenotype info", td(pdata()))
      incProgress(1/5, detail = "writing geneset annotation")
      writeData(wb, sheet = "Geneset annot", attr(fdata(), 
                                                  "GS"))
      incProgress(1/5, detail = "Saving table")
      saveWorkbook(wb, file = file, overwrite = TRUE)
    })
  })
  output$summary <- renderUI({
    if (!vEset()) {
      txt <- sprintf("<h1 style=\"display:inline;\">%s</h1> <h3 style=\"display:inline;\"><sup>%s</sup></h3>", 
                     appName, paste0("v", appVersion))
    }
    else {
      txt <- sprintf("<h1 style=\"display:inline;\">%s</h1> <h3 style=\"display:inline;\"><sup>%s</sup>  --   %s features and %s samples:</h3>", 
                     appName, paste0("v", appVersion), nrow(expr()), 
                     ncol(expr()))
    }
    HTML(txt)
  })
  v1 <- callModule(L1_data_space_module, id = "dataspace", 
                   expr = expr, pdata = pdata, fdata = fdata, reactive_x_s = d_s_x, 
                   reactive_y_s = d_s_y, reactive_x_f = d_f_x, reactive_y_f = d_f_y, 
                   status = esv_status, cormat = cormat , imputationMethod = imputationMethod)
  sameValues <- function(a, b) {
    if (is.null(a) || is.null(b)) 
      return(FALSE)
    all(sort(a) == sort(b))
  }
  ri <- reactiveVal()
  observeEvent(v1(), {
    ri(c(v1()$feature))
  })
  observeEvent(expr(), ri(NULL))
  rh <- reactiveVal()
  observeEvent(v1(), {
    rh(c(v1()$sample))
  })
  observeEvent(expr(), rh(NULL))
  v2 <- callModule(L1_result_space_module, id = "resultspace", 
                   reactive_expr = expr, reactive_phenoData = pdata, reactive_featureData = fdata, 
                   reactive_i = ri, reactive_highlight = rh, additionalTabs = additionalTabs, 
                   object = reactive_eset, status = esv_status)
  dir <- reactiveVal()
  observe({
    dd <- getwd()
    if (!is.null(.dir())) 
      dd <- .dir()
    dir(dd)
  })
  savedSS <- reactiveVal()
  observe({
    req(.dir())
    if (is.null(input$selectFile) || nchar(input$selectFile) == 
        0) 
      fs <- "ESVObj.RDS"
    else fs <- input$selectFile
    fl <- paste0("ESVSnapshot_", fs, "_")
    ff <- list.files(.dir(), pattern = fl)
    if (length(ff) == 0) 
      return(NULL)
    r <- sub(fl, "", ff)
    r <- sub(".ESS$", "", r)
    df <- data.frame(name = r, link = ff, stringsAsFactors = FALSE, 
                     check.names = FALSE)
    savedSS(df)
  })
  shinyInput <- function(FUN, len, id, ...) {
    inputs <- c()
    for (i in len) {
      inputs <- c(inputs, as.character(FUN(paste0(id, i), 
                                           ...)))
    }
    inputs
  }
  output$tab_saveSS <- renderDT({
    req(nrow(dt <- savedSS()) > 0)
    dt$delete <- shinyInput(actionButton, dt$name, "deletess_", 
                            label = "Delete", onclick = sprintf("Shiny.setInputValue(\"%s\",  this.id)", 
                                                                ns("deletess_button")))
    DT::datatable(dt[, c(1, 3), drop = FALSE], rownames = FALSE, 
                  colnames = c(NULL, NULL, NULL), selection = list(mode = "single", 
                    target = "cell", selectable = -cbind(seq_len(nrow(dt)), 
                    1)), escape = FALSE, options = list(dom = "t", 
                    autoWidth = FALSE, style = "compact-hover", scrollY = "450px", 
                    paging = FALSE, columns = list(list(width = "85%"), 
                                                    list(width = "15%"))))
  })
  selectedSS <- reactiveVal()
  observe({
    ss <- input$tab_saveSS_cells_selected

    if (length(ss) == 0 || ss[2] > 0) 
      return(NULL)
    selectedSS(ss[1])
  })
  observeEvent(list(v1(), v2()), {
    selectedSS(NULL)
  })
  deleteSS <- reactiveVal()
  observeEvent(input$deletess_button, {
    selectedRow <- sub("deletess_", "", input$deletess_button)
    deleteSS(selectedRow)
  })
  observeEvent(input$snapshot, {
    showModal(modalDialog(title = NULL, fluidRow(column(9, 
                                                        textInput(ns("snapshot_name"), label = "Save new snapshot", 
                                                                  placeholder = "snapshot name", width = "100%")), 
                                                 column(3, style = "padding-top:25px", actionButton(ns("snapshot_save"), 
                                                                                                    label = "Save")), ), hr(), strong("Load saved snapshots:"), 
                          DTOutput(ns("tab_saveSS")), footer = NULL, easyClose = TRUE))
  })
  observeEvent(input$snapshot_save, {
    if (is.null(input$selectFile) || nchar(input$selectFile) == 
        0) 
      fs <- "ESVObj.RDS"
    else fs <- input$selectFile
    df <- savedSS()
    if (!is.null(df)) {
      if (input$snapshot_name %in% df$name) {
        showModal(modalDialog(title = "FAILED!", "Snapshot with this name already exists, please give a different name."))
        return(NULL)
      }
    }
    obj <- c(attr(v1(), "status"), v2(), active_feature = list(ri()), 
             active_sample = list(rh()))
    flink <- file.path(.dir(), paste0("ESVSnapshot_", fs, 
                                      "_", input$snapshot_name, ".ESS"))
    saveRDS(obj, flink)
    df <- rbind(df, data.frame(name = input$snapshot_name, 
                               link = basename(flink)), stringsAsFactors = FALSE)
    dt <- df[order(df$name), ]
    savedSS(dt)
    removeModal()
  })
  esv_status <- reactiveVal()
  observeEvent(selectedSS(), {
    req(nrow(df <- savedSS()) > 0)
    if (length(i <- selectedSS()) == 0) 
      return(NULL)
    removeModal()
    esv_status(NULL)
    esv_status(readRDS(file.path(.dir(), df[i, 2])))
  })
  observeEvent(deleteSS(), {
    req(nrow(df <- savedSS()) > 0)
    req(i <- match(deleteSS(), df$name))
    df <- savedSS()
    unlink(file.path(.dir(), df[i, 2]))
    df <- df[-i, , drop = FALSE]
    savedSS(df)
  })
}


app_ui <- function (id, showDropList = TRUE, activeTab = "Feature") 
{
  ns <- NS(id)
  comp <- list(useShinyjs(), style = "background:white;", absolutePanel(top = 5, 
                                                                        right = 20, style = "z-index: 9999;", width = 115, downloadButton(outputId = ns("download"), 
                                                                                label = "xlsx", class = NULL), actionButton(ns("snapshot"), 
                                                                                label = NULL, icon = icon("camera-retro"))), shinyjs::hidden(div(id = ns("contents"), 
                                                                                column(6, L1_data_space_ui(ns("dataspace"), activeTab = activeTab)), 
                                                                                column(6, L1_result_space_ui(ns("resultspace"))))))
  if (showDropList) {
    l2 <- list(shinycssloaders::withSpinner(uiOutput(ns("summary")), 
                                            hide.ui = FALSE, type = 8, color = "green"), br(), 
               absolutePanel(top = 8, right = 140, style = "z-index: 9999;", 
                             selectizeInput(inputId = ns("selectFile"), label = NULL, 
                                            choices = NULL, width = "500px", options = list(placeholder = "Select a dataset here"))))
    comp <- c(l2, comp)
  }
  do.call(fluidRow, comp)
}


asEsetWithAttr <- function (x) 
{
  if (inherits(x, "SummarizedExperiment")) {
    eset <- as(x, "ExpressionSet")
    colnames(pData(eset)) <- colnames(colData(x))
    colnames(fData(eset)) <- colnames(rowData(x))
    DFattrs <- c("rownames", "nrows", "listData", "elementType", 
                 "elementMetadata", "metadata", "class")
    for (i in setdiff(names(attributes(colData(x))), DFattrs)) attr(pData(eset), 
                                                                    i) <- attr(colData(x), i)
    for (i in setdiff(names(attributes(rowData(x))), DFattrs)) attr(fData(eset), 
                                                                    i) <- attr(rowData(x), i)
    SEattrs <- c("assays", "colData", "NAMES", "elementMetadata", 
                 "metadata", "class")
    for (i in setdiff(names(attributes(x)), SEattrs)) attr(eset, 
                                                           i) <- attr(x, i)
  }
  else if (inherits(x, "ExpressionSet")) {
    eset <- x
  }
  else stop("x should be either an SummarizedExperiment or ExpressionSet")
  eset
}


attr4selector_module <- function (input, output, session, reactive_meta = reactive(NULL), 
                                  reactive_expr = reactive(NULL), reactive_triset = reactive(NULL), 
                                  pre_volcano = reactive(FALSE), reactive_status = reactive(NULL)) 
{
  ns <- session$ns
  params <- reactiveValues(highlight = NULL, highlightName = NULL, 
                           color = NULL, shape = NULL, size = NULL, tooltips = NULL, 
                           cutoff = NULL)
  selectColor_s1 <- reactiveVal()
  selectColor_s2 <- reactiveVal()
  selectColor_s3 <- reactiveVal()
  selectShape_s1 <- reactiveVal()
  selectShape_s2 <- reactiveVal()
  selectShape_s3 <- reactiveVal()
  selectSize_s1 <- reactiveVal()
  selectSize_s2 <- reactiveVal()
  selectSize_s3 <- reactiveVal()
  selectTooltip_s1 <- reactiveVal()
  selectTooltip_s2 <- reactiveVal()
  selectTooltip_s3 <- reactiveVal()
  searchOnCol_s1 <- reactiveVal()
  searchOnCol_s2 <- reactiveVal()
  searchOnCol_s3 <- reactiveVal()
  selectColor <- callModule(triselector_module, id = "selectColorUI", 
                            reactive_x = reactive_triset, label = "Color", reactive_selector1 = selectColor_s1, 
                            reactive_selector2 = selectColor_s2, reactive_selector3 = selectColor_s3)
  selectShape <- callModule(triselector_module, id = "selectShapeUI", 
                            reactive_x = reactive_triset, label = "Shape", reactive_selector1 = selectShape_s1, 
                            reactive_selector2 = selectShape_s2, reactive_selector3 = selectShape_s3)
  selectSize <- callModule(triselector_module, id = "selectSizeUI", 
                           reactive_x = reactive_triset, label = "Size", reactive_selector1 = selectSize_s1, 
                           reactive_selector2 = selectSize_s2, reactive_selector3 = selectSize_s3)
  selectTooltip <- callModule(triselector_module, id = "selectTooltipUI", 
                              reactive_x = reactive_triset, label = "Tooltips", reactive_selector1 = selectTooltip_s1, 
                              reactive_selector2 = selectTooltip_s2, reactive_selector3 = selectTooltip_s3)
  searchOnCol <- callModule(triselector_module, id = "selectSearchCol", 
                            reactive_x = reactive_triset, label = "Search", reactive_selector1 = searchOnCol_s1, 
                            reactive_selector2 = searchOnCol_s2, reactive_selector3 = searchOnCol_s3)
  vv <- reactive(varSelector(searchOnCol(), expr = reactive_expr(), 
                             meta = reactive_meta()))
  pre_search <- reactiveVal()
  observe({
    updateSelectInput(session, "searchon", choices = vv(), 
                      selected = pre_search())
  })
  observe(updateCheckboxInput(session, "showSearchBox", value = !is.null(vv())))
  val_xcut <- reactive({
    text2num(input$xcut)
  })
  val_ycut <- reactive({
    text2num(input$ycut)
  })
  observeEvent(list(val_xcut(), val_ycut()), {
    if (is.numeric(val_xcut()) && is.null(val_ycut())) {
      ac <- c("None", "left", "right")
    }
    else if (is.null(val_xcut()) && is.numeric(val_ycut())) {
      ac <- c("None", "top", "bottom")
    }
    else if (is.numeric(val_xcut()) && is.numeric(val_ycut())) {
      ac <- c("None", "volcano", "left", "right", "top", 
              "bottom", "topleft", "topright", "bottomleft", 
              "bottomright")
    }
    else ac <- "None"
    ps <- ac[1]
    if (!is.null(input$scorner) && input$scorner %in% ac) 
      ps <- input$scorner
    updateSelectInput(session, inputId = "scorner", choices = ac, 
                      selected = ps)
  })
  observeEvent(pre_volcano(), {
    if (pre_volcano()) {
      l <- list(x = val_xcut(), y = val_ycut(), corner = "volcano")
      attr(l, "seed") <- Sys.time()
      params$cutoff <- l
      updateSelectInput(session, inputId = "scorner", selected = "volcano")
    }
    else {
      params$cutoff <- list(x = val_xcut(), y = val_ycut(), 
                            corner = "None")
      updateSelectInput(session, inputId = "scorner", selected = "None")
    }
  })
  searchValue <- reactiveVal()
  observe({
    foo <- function() searchValue(input$searchon)
    debounce(foo, 1000)
  })
  observe({
    if (is.null(vv())) {
      updateSelectInput(session, "searchon", choices = NULL, 
                        selected = NULL)
      searchValue(NULL)
      pre_search(NULL)
    }
  })
  observe({
    req(vv())
    req(searchValue())
    params$highlight <- which(vv() %in% searchValue())
    isolate(params$highlightName <- searchOnCol()$variable)
  })
  observe(params$color <- varSelector(selectColor(), reactive_expr(), 
                                      reactive_meta(), alternative = selectShape()$variable))
  observe(params$shape <- varSelector(selectShape(), reactive_expr(), 
                                      reactive_meta(), alternative = selectColor()$variable))
  observe(params$size <- varSelector(selectSize(), reactive_expr(), 
                                     reactive_meta()))
  observe(params$tooltips <- varSelector(selectTooltip(), reactive_expr(), 
                                         reactive_meta()))
  acorner <- reactiveVal()
  i_xcut <- reactiveVal()
  i_ycut <- reactiveVal()
  observe({
    req(!is.null(input$scorner) && nchar(input$scorner) != 
          0)
    acorner(input$scorner)
    i_xcut(input$xcut)
    i_ycut(input$ycut)
    l <- list(x = val_xcut(), y = val_ycut(), corner = input$scorner)
    attr(l, "seed") <- Sys.time()
    params$cutoff <- l
  })
  observe({
    params$status <- list(selectColor = selectColor(), selectShape = selectShape(), 
                          selectSize = selectSize(), selectTooltip = selectTooltip(), 
                          searchOnCol = searchOnCol(), searchValue = searchValue(), 
                          xcut = i_xcut(), ycut = i_ycut(), acorner = acorner())
  })
  observe({
    if (is.null(s <- reactive_status())) 
      return(NULL)
    selectColor_s1(s$selectColor[[1]])
    selectColor_s2(s$selectColor[[2]])
    selectColor_s3(s$selectColor[[3]])
  })
  observe({
    if (is.null(s <- reactive_status())) 
      return(NULL)
    selectShape_s1(s$selectShape[[1]])
    selectShape_s2(s$selectShape[[2]])
    selectShape_s3(s$selectShape[[3]])
  })
  observe({
    if (is.null(s <- reactive_status())) 
      return(NULL)
    selectSize_s1(s$selectSize[[1]])
    selectSize_s2(s$selectSize[[2]])
    selectSize_s3(s$selectSize[[3]])
  })
  observe({
    if (is.null(s <- reactive_status())) 
      return(NULL)
    selectTooltip_s1(s$selectTooltip[[1]])
    selectTooltip_s2(s$selectTooltip[[2]])
    selectTooltip_s3(s$selectTooltip[[3]])
  })
  observe({
    if (is.null(s <- reactive_status())) 
      return(NULL)
    searchOnCol_s1(s$searchOnCol[[1]])
    searchOnCol_s2(s$searchOnCol[[2]])
    searchOnCol_s3(s$searchOnCol[[3]])
  })
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status())) 
      return(NULL)
    updateTextInputIcon(session, "xcut", value = s$xcut)
    updateTextInputIcon(session, "ycut", value = s$ycut)
    updateSelectInput(session, "scorner", selected = s$acorner)
  })
  observe({
    if (is.null(s <- reactive_status())) 
      return(NULL)
    pre_search(s$searchValue)
  })
  params
}


attr4selector_ui <- function (id, circle = TRUE, right = FALSE) 
{
  ns <- NS(id)
  dropdown(margin = "25px", circle = circle, right = right, 
           status = "default", icon = icon("cog"), width = "788px", 
           tooltip = tooltipOptions(title = "Click to modify figure!"), 
           br(), triselector_ui(ns("selectColorUI")), triselector_ui(ns("selectShapeUI")), 
           triselector_ui(ns("selectSizeUI")), triselector_ui(ns("selectTooltipUI")), 
           triselector_ui(ns("selectSearchCol")), conditionalPanel("1 == 2", 
          checkboxInput(ns("showSearchBox"), label = "show", 
          value = FALSE)), conditionalPanel("input.showSearchBox == true", 
          ns = ns, div(style = "padding-left:100px; padding-right:0px; padding-top:0px; padding-bottom:0px", 
          selectInput(ns("searchon"), label = NULL, choices = NULL, 
          multiple = TRUE, width = "100%"))), fluidRow(column(2), 
          column(4, offset = 0, style = "padding:2px;", textInputIcon(inputId = ns("xcut"), 
          label = "Select points by x/y cutoffs", value = "log10(2)", 
          placeholder = "e.g. -1 or -log10(2)", icon = list("x-cut"))), 
          column(4, offset = 0, style = "padding-left:5px; padding-right:2px; padding-top:27px; padding-bottom:2px;", 
          textInputIcon(inputId = ns("ycut"), label = NULL, 
          value = "-log10(0.05)", placeholder = "e.g 2 or -log10(0.05)", 
          icon = list("y-cut"))), column(2, offset = 0, 
          style = "padding-left:5px; padding-right:20px; padding-top:4px; padding-bottom:2px;", 
          selectInput(inputId = ns("scorner"), label = "Area", 
              choices = "None", selectize = TRUE))))

}


correlationAnalysis <- function (x, pheno, min.value = 12, prefix = "Cor") 
{
  if (is.data.frame(pheno)) {
    cn <- colnames(pheno)
    rn <- rownames(pheno)
    pheno <- apply(pheno, 2, function(x) as.numeric(as.character(x)))
    colnames(pheno) <- cn
    rownames(pheno) <- rn
  }
  if (is.null(colnames(pheno))) 
    colnames(pheno) <- make.names(seq_len(ncol(pheno)))
  if (ncol(x) != nrow(pheno)) 
    stop("ncol x != nrow pheno")
  r <- psych::corr.test(t(x), pheno, use = "pair", adjust = "none", 
                        ci = FALSE)
  rs <- lapply(seq_len(ncol(pheno)), function(i) {
    if (length(r$n) == 1) 
      nv <- r$n
    else nv <- r$n[, i]
    df <- data.frame(R = r$r[, i], N = nv, P = r$p[, i], 
                     logP = -log10(r$p[, i]))
    texp <- x[, !is.na(pheno[, i])]
    df$range <- rowMaxs(texp, na.rm = TRUE) - rowMins(texp, 
                                                      na.rm = TRUE)
    ii <- which(df$N < min.value)
    df[c("R", "P", "logP")] <- lapply(df[c("R", "P", "logP")], 
                                      function(x) {
                                        x[ii] <- NA
                                        x
                                      })
    if (all(is.na(df$R))) 
      return(df[, character(0)])
    colnames(df) <- paste(colnames(pheno)[i], colnames(df), 
                          sep = "|")
    df
  })
  rs <- do.call(cbind, rs)
  if (ncol(rs) > 0) 
    colnames(rs) <- paste(prefix, colnames(rs), sep = "|")
  rs
}


csc2list <- function (x) 
{
  if (is.matrix(x)) {
    if (is.null(rownames(x)) || is.null(colnames(x))) 
      stop("csmlist: x need to have dimnames!")
    x[x == 0] <- NA
    df <- reshape2::melt(x, na.rm = TRUE)
    colnames(df) <- c("featureId", "gsId", "weight")
  }
  else if (inherits(x, "dgCMatrix")) {
    if (is.null(x@Dimnames[[1]]) || is.null(x@Dimnames[[2]])) 
      stop("csmlist: x need to have dimnames!")
    df <- data.frame(featureId = x@Dimnames[[1]][x@i + 1], 
                     gsId = rep(x@Dimnames[[2]], x@p[-1] - x@p[-length(x@p)]), 
                     stringsAsFactors = TRUE)
    if (hasAttr(x, "x")) 
      df$weight <- x@x
  }
  else stop("x should be either a matrix or CsparseMatrix")
  rownames(df) <- NULL
  df
}


dataTable_module <- function (input, output, session, reactive_data, selector = TRUE, 
                              columns = NULL, tab_status = reactive(NULL), tab_rows = reactive(NULL)) 
{
  ns <- session$ns
  selectedRowOrCol <- reactiveVal(TRUE)
  notNullAndPosLength <- function(x) !is.null(x) && length(x) > 
    0
  observeEvent(reactive_data(), selectedRowOrCol(TRUE))
  observeEvent(input$clear, selectedRowOrCol(TRUE))
  observe(selectedRowOrCol(tab_rows()))
  rdd <- reactive({
    req(reactive_data())
    if (is.matrix(reactive_data())) {
      x <- as.data.frame(reactive_data())
    }
    else if (is.data.frame(reactive_data())) 
      x <- reactive_data()
    else stop("reactive_data shold be either a matrix or data.frame")
    x[selectedRowOrCol(), , drop = FALSE]
  })
  callModule(dataTableDownload_module, id = "downloadTable", 
             reactive_table = rdd, prefix = "viewerTable_")
  cols <- eventReactive(reactive_data(), {
    cn <- intersect(columns, colnames(rdd()))
    if (length(cn) == 0) 
      cn <- grep("^General\\|", colnames(rdd()), ignore.case = TRUE, 
                 value = TRUE)
    if (length(cn) == 0) 
      cn <- colnames(rdd())
    opt <- NULL
    if (selector) {
      optx <- setdiff(colnames(rdd()), cn)
      if (length(optx) > 0) 
        opt <- str_split_fixed(optx, pattern = "\\|", 
                               n = 3)
    }
    list(shown = cn, opt = opt)
  })
  addcols <- callModule(triselector_module, id = "select", 
                        reactive_x = reactive({
                          req(cols()$opt)
                          req(nrow(cols()$opt) > 0)
                          cols()$opt
                        }), label = "Add column")
  scn <- reactiveVal(NULL)
  observe(scn(cols()$shown))
  observeEvent(addcols(), {
    req(!addcols()$variable %in% c("", "Select a variable!"))
    oc <- scn()
    nc <- unique(c(oc, paste(addcols(), collapse = "|")))
    nc <- intersect(nc, colnames(rdd()))
    req(nc)
    scn(nc)
  })
  output$selector <- renderUI({
    req(cols()$opt)
    req(nrow(cols()$opt) > 0)
    triselector_ui(ns("select"))
  })
  formatTab <- function(tab, sel) {
    ci <- unname(which(vapply(tab, inherits, c("factor", 
                                               "character"), FUN.VALUE = logical(1))))
    if (length(ci) > 0) 
      tab[ci] <- lapply(tab[ci], function(x) {
        x[is.na(x)] <- ""
        x
      })
    dt <- DT::datatable(tab, selection = list(mode = c("single", 
                                                       "multiple")[as.integer(sel) + 1], selected = tab_status()$rows_selected, 
                                              target = "row"), rownames = FALSE, filter = "top", 
                        class = "table-bordered compact nowrap", options = list(scrollX = TRUE, 
                                                                                pageLength = 25, dom = "tip", columnDefs = list(list(targets = ci - 
                                                                                                                                       1, render = DT::JS("function(data, type, row, meta) {", 
                                                                                                                                                          "return type === 'display' && data.length > 50 ?", 
                                                                                                                                                          "'<span title=\"' + data + '\">' + data.substr(0, 50) + '...</span>' : data;", 
                                                                                                                                                          "}"))), stateSave = TRUE, stateDuration = -1, 
                                                                                searchCols = getSearchCols(tab_status()), order = getOrderCols(tab_status()), 
                                                                                displayStart = tab_status()$start))
    DT::formatStyle(dt, columns = seq_len(ncol(tab)), fontSize = "90%")
  }
  observeEvent(tab_status(), {
    if (!is.null(i <- tab_status()$showColumns)) 
      scn(i)
    updateSwitchInput(session, "multisel", value = tab_status()$multiSelection)
  })
  output$table <- DT::renderDataTable({
    req(scn())
    tab <- rdd()[, scn(), drop = FALSE]
    i <- which(vapply(tab, function(x) is.numeric(x) && !is.integer(x), 
                      logical(1)))
    if (any(i)) 
      tab[i] <- lapply(tab[i], round, digits = 4)
    formatTab(tab, sel = input$multisel)
  })
  tabproxy <- dataTableProxy(ns("table"))
  eventReactive(list(input$table_rows_selected, input$table_state), 
                {
                  r <- character(0)
                  if (!is.null(tab_rows())) 
                    r <- tab_rows()
                  if (notNullAndPosLength(input$table_rows_selected)) 
                    r <- rownames(rdd())[input$table_rows_selected]
                  sta <- input$table_state
                  sta$showColumns <- scn()
                  sta$multiSelection <- input$multisel
                  sta$rows_selected <- input$table_rows_selected
                  attr(r, "status") <- sta
                  r
                })
}


dataTable_ui <- function (id) 
{
  ns <- NS(id)
  tagList(fluidRow(column(3, actionButton(ns("clear"), "Show all")), 
                   column(6, align = "center", shinyWidgets::switchInput(inputId = ns("multisel"), 
                                                                         label = "Multiple_selection", labelWidth = "125px")), 
                   column(3, dataTableDownload_ui(ns("downloadTable"), showTable = FALSE), 
                          align = "right")), uiOutput(ns("selector")), DT::dataTableOutput(ns("table")))
}


dataTableDownload_module <- function (input, output, session, reactive_table, tab_status = NULL, 
                                      reactive_cols = reactive(NULL), prefix = "", pageLength = 10, 
                                      sortBy = NULL, decreasing = TRUE) 
{
  ns <- session$ns
  rtab <- reactive({
    req(tt <- reactive_table())
    if (is.matrix(tt)) 
      tt <- as.data.frame(tt, stringsAsFactors = FALSE)
    tt
  })
  output$downloadData <- downloadHandler(filename = function() {
    paste0(prefix, Sys.time(), ".tsv")
  }, content = function(file) {
    tab <- rtab()
    ic <- which(vapply(tab, is.list, logical(1)))
    if (length(ic) > 0) {
      for (ii in ic) {
        tab[, ii] <- vapply(tab[, ii], paste, collapse = ";", 
                            FUN.VALUE = character(1))
      }
    }
    write.table(tab, file, col.names = TRUE, row.names = FALSE, 
                quote = FALSE, sep = "\t")
  })
  output$showButton <- renderUI({
    req(rtab())
    downloadLink(ns("downloadData"), "Save table")
  })
  formatTab <- function(tab, sel = 0, pageLength = 10) {
    dt <- DT::datatable(tab, selection = c("single", "multiple")[as.integer(sel) + 
                                                                   1], rownames = FALSE, filter = "top", class = "table-bordered compact nowrap", 
                        options = list(scrollX = TRUE, pageLength = pageLength, 
                                       dom = "tip", stateSave = TRUE, stateDuration = -1, 
                                       searchCols = getSearchCols(tab_status), order = getOrderCols(tab_status)))
    DT::formatStyle(dt, columns = seq_len(ncol(tab)), fontSize = "90%")
  }
  tabsort <- reactive({
    req(tab <- rtab())
    index <- seq_len(nrow(tab))
    if (!is.null(sortBy)) {
      if (sortBy %in% colnames(tab)) {
        o <- order(tab[, sortBy], decreasing = decreasing)
        tab <- tab[o, ]
        index <- index[o]
      }
    }
    if (!is.null(reactive_cols())) 
      tab <- tab[, reactive_cols()]
    ic <- which(vapply(tab, function(x) is.numeric(x) & !is.integer(x), 
                       logical(1)))
    if (length(ic) > 0) 
      tab[, ic] <- lapply(tab[, ic, drop = FALSE], signif, 
                          digits = 3)
    list(tab = tab, index = index)
  })
  output$table <- DT::renderDataTable(formatTab(tabsort()$tab, 
                                                pageLength = pageLength))
  reactive({
    ii <- tabsort()$index[input$table_rows_selected]
    attr(ii, "status") <- input$table_state
    ii
  })
}


dataTableDownload_ui <- function (id, showTable = TRUE) 
{
  ns <- NS(id)
  if (showTable) {
    r <- tagList(uiOutput(ns("showButton")), DT::dataTableOutput(ns("table")))
  }
  else r <- uiOutput(ns("showButton"))
  r
}


downloadUPRefProteome <- function (id, domain = c("Eukaryota", "Archaea", "Bacteria", 
                                                  "Viruses")[1], destdir = "./") 
{
  url <- sprintf("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/%s/%s", 
                 domain, id)
  url2 <- curl(url)
  on.exit(close(url2))
  con <- readLines(url2)
  con <- na.omit(stringr::str_match(con, paste0(id, "_s*(.*?)\\s*.dat.gz"))[, 
                                                                            1])
  con <- con[!grepl("additional", con)]
  durl <- paste(url, con, sep = "/")
  download.file(durl, destfile = file.path(destdir, con))
  con
}


draw_roc_pr <- function (value, label) 
{
  pred_one <- function(value, label) {
    cal_auc <- function(pred) {
      perf <- performance(pred, measure = "auc")
      perf@y.values[[1]]
    }
    pred <- prediction(value, label)
    auc <- cal_auc(pred)
    if (auc < 0.5) {
      label <- 1 - label
      pred <- prediction(value, label)
      auc <- cal_auc(pred)
    }
    list(roc = performance(pred, measure = "tpr", x.measure = "fpr"), 
         auc = auc, pr = performance(pred, measure = "prec", 
                                     x.measure = "rec"))
  }
  label <- as.character(label)
  i <- !is.na(value) & !is.na(label)
  value <- value[i]
  label <- label[i]
  nlab <- length(unique(label))
  if (nlab == 1 || nlab > 12) {
    message("label should have 2 - 12 unique values!")
    plot(0, col = NA, xlab = "", ylab = "", axes = FALSE)
    mtext(text = "label should have 2 - 12 unique values!", 
          side = 3, line = -5)
    return()
  }
  dum <- model.matrix(~label - 1)
  if (ncol(dum) == 2) {
    dum <- dum[, 1, drop = FALSE]
    colnames(dum) <- "Label"
  }
  else colnames(dum) <- sub("label", "", colnames(dum))
  curves <- lapply(seq_len(ncol(dum)), function(i) pred_one(value, 
                                                            dum[, i]))
  layout(matrix(1:2, 1, 2))
  cols <- nColors(ncol(dum))
  for (i in seq_along(curves)) plot(curves[[i]]$roc, add = i > 
                                      1, col = cols[i], main = "ROC Curve", lwd = 2)
  lab <- paste0(colnames(dum), " | ", colSums(dum), " | ", 
                signif(sapply(curves, "[[", "auc"), 2))
  legend("bottomright", col = cols, legend = lab, lty = 1, 
         lwd = 2, bty = "n", title = "label | n | AUC")
  for (i in seq_along(curves)) plot(curves[[i]]$pr, add = i > 
                                      1, col = cols[i], xlim = c(0, 1), ylim = c(0, 1), main = "PR Curve", 
                                    lwd = 2)
}


drawButton <- function (id) 
{
  list(name = "Draw", icon = list(path = "M109.46 244.04l134.58-134.56-44.12-44.12-61.68 61.68a7.919 7.919 0 0 1-11.21 0l-11.21-11.21c-3.1-3.1-3.1-8.12 0-11.21l61.68-61.68-33.64-33.65C131.47-3.1 111.39-3.1 99 9.29L9.29 99c-12.38 12.39-12.39 32.47 0 44.86l100.17 100.18zm388.47-116.8c18.76-18.76 18.75-49.17 0-67.93l-45.25-45.25c-18.76-18.76-49.18-18.76-67.95 0l-46.02 46.01 113.2 113.2 46.02-46.03zM316.08 82.71l-297 296.96L.32 487.11c-2.53 14.49 10.09 27.11 24.59 24.56l107.45-18.84L429.28 195.9 316.08 82.71zm186.63 285.43l-33.64-33.64-61.68 61.68c-3.1 3.1-8.12 3.1-11.21 0l-11.21-11.21c-3.09-3.1-3.09-8.12 0-11.21l61.68-61.68-44.14-44.14L267.93 402.5l100.21 100.2c12.39 12.39 32.47 12.39 44.86 0l89.71-89.7c12.39-12.39 12.39-32.47 0-44.86z", 
                                  width = 1000, height = 1000, transform = "matrix(0.75 0 0 -0.75 0 1000) scale(2.5)"), 
       click = htmlwidgets::JS(paste0("function(gd) {\n             // priority 'event' causes invalidation irrespective of value\n             Shiny.setInputValue('", 
                                      id, "', 'bar', {priority: 'event'});\n          }")))
}


enrichment_analysis_module <- function (input, output, session, reactive_featureData, reactive_i) 
{
  ns <- session$ns
  reactive_pathway <- reactive({
    attr(reactive_featureData(), "GS")
  })
  rii <- reactiveVal()
  observe({
    req(reactive_i())
    if (length(reactive_i()) <= 1) 
      return(NULL)
    if (length(reactive_i()) <= 3) 
      rii("notest")
    else rii(reactive_i())
  })
  oraTab <- reactive({
    req(rii())
    notest <- "No geneset has been tested, please try to include more input feature IDs!"
    if (rii()[1] == "notest") 
      return(notest)
    tab <- vectORATall(reactive_pathway(), i = rii(), background = nrow(reactive_featureData()))
    if (is.null(tab)) 
      return(notest)
    ic <- which(vapply(tab, function(x) is.numeric(x) & !is.integer(x), 
                       logical(1)))
    tab[, ic] <- lapply(tab[, ic], signif, digits = 3)
    tab <- tab[which(tab$p.adjusted < 0.1 | tab$p.value < 
                       0.05 | tab$OR >= 3), ]
    tab
  })
  output$errorMsg <- renderText({
    req(is.character(oraTab()))
    oraTab()
  })
  output$error <- renderUI(verbatimTextOutput(ns("errorMsg")))
  vi <- callModule(dataTableDownload_module, id = "stab", reactive_table = reactive({
    req(is.data.frame(oraTab()))
    oraTab()
  }), reactive_cols = reactive(setdiff(colnames(oraTab()), 
                                       "overlap_ids")), prefix = "ORA_", sortBy = "p.value", 
  decreasing = FALSE)
  hd <- reactive({
    req(is.data.frame(oraTab()))
    req(i <- vi())
    ii <- grep("^General", colnames(reactive_featureData()), 
               ignore.case = TRUE)
    if (length(ii) == 0) 
      ii <- seq_len(min(3, ncol(reactive_featureData())))
    i <- oraTab()[i, ]
    hid <- i$overlap_ids[[1]]
    req(hid)
    df1 <- reactive_featureData()[hid, ii, drop = FALSE]
    df1 <- cbind(Overlap = "+", df1)
    apath <- reactive_pathway()[reactive_pathway()$gsId == 
                                  i$pathway, ]
    aid <- setdiff(apath$featureId, hid)
    if (length(aid) > 0) {
      df2 <- reactive_featureData()[aid, ii, drop = FALSE]
      df2 <- cbind(Overlap = "", df2)
      df1 <- rbind(df1, df2)
    }
    df1
  })
  vi2 <- callModule(dataTableDownload_module, id = "overlapTab", 
                    reactive_table = hd, prefix = "ORA_overlapGenes_")
}


enrichment_analysis_ui <- function (id) 
{
  ns <- NS(id)
  tagList(uiOutput(ns("error")), dataTableDownload_ui(ns("stab")), 
          dataTableDownload_ui(ns("overlapTab")))
}


enrichment_fgsea_module <- function (input, output, session, reactive_featureData, reactive_status = reactive(NULL)) 
{
  ns <- session$ns
  triset <- reactive({
    fd <- reactive_featureData()
    cn <- colnames(fd)[vapply(fd, is.numeric, logical(1)) & 
                         !grepl("^GS\\|", colnames(fd))]
    str_split_fixed(cn, "\\|", n = 3)
  })
  xax <- reactiveVal()
  v1 <- callModule(triselector_module, id = "tris_fgsea", reactive_x = triset, 
                   label = "Input variable", reactive_selector1 = reactive(xax()$v1), 
                   reactive_selector2 = reactive(xax()$v2), reactive_selector3 = reactive(xax()$v3))
  gsInfo <- reactive({
    fdgs <- attr(reactive_featureData(), "GS")
    uniqueGs <- unique(fdgs$gsId)
    names(uniqueGs) <- uniqueGs
    list(gs = fdgs, desc = uniqueGs)
  })
  tab <- reactive({
    req(!v1()$variable %in% c("Select a variable!", ""))
    scc <- paste(v1(), collapse = "|")
    req(scc %in% colnames(reactive_featureData()))
    stats <- reactive_featureData()[, scc]
    names(stats) <- rownames(reactive_featureData())
    stats <- na.omit(stats)
    fdgs <- gsInfo()$gs[gsInfo()$gs$featureId %fin% names(stats), 
    ]
    if (nrow(fdgs) < 3) {
      message("Perhaps a problem ... enrichment_fgsea_module")
      return(NULL)
    }
    res <- fgsea1(fdgs, stats = stats, minSize = 3, maxSize = 500, 
                  gs_desc = gsInfo()$desc)
    cn <- colnames(res)
    cn[cn == "ES"] <- "enrichment score (ES)"
    cn[cn == "NES"] <- "normalized ES"
    cn[cn == "desc"] <- "description"
    colnames(res) <- cn
    list(pathway_mat = fdgs, table = res[order(abs(res$"normalized ES"), 
                                               decreasing = TRUE), , drop = FALSE], stats = stats, 
         statsNames = names(stats))
  })
  vi <- callModule(dataTableDownload_module, id = "stab", reactive_table = reactive(tab()$table), 
                   reactive_cols = reactive(setdiff(colnames(tab()$table), 
                                                    "leadingEdge")), prefix = "fgsea_")
  output$bplot <- renderPlotly({
    hid <- bid <- NULL
    if (!is.null(i <- vi()) && length(vi()) > 0) {
      i <- tab()$table[i, ]
      hid <- i$leadingEdge[[1]]
      bid <- setdiff(tab()$pathway_mat$featureId[tab()$pathway_mat$gsId == 
                                                   i$pathway], hid)
      if (length(bid) == 0) 
        bid <- NULL
      hid <- fmatch(hid, tab()$statsNames)
      bid <- fmatch(bid, tab()$statsNames)
    }
    plotly_barplot(x = tab()$stats, names = tab()$statsNames, 
                   highlight = hid, highlight_color = "red", highlight_width = 2, 
                   highlight_legend = "Leading edges", background = bid, 
                   background_color = "gray", background_width = 2, 
                   background_legend = "background", ylab = "Rankding stats", 
                   xlab = "", sort = "dec", source = ns("plotlybarchart"))
  })
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status())) 
      return()
    xax(NULL)
    xax(list(v1 = s$xax[[1]], v2 = s$xax[[2]], v3 = s$xax[[3]]))
  })
  reactive(list(xax = v1()))
}


enrichment_fgsea_ui <- function (id) 
{
  ns <- NS(id)
  tagList(triselector_ui(ns("tris_fgsea"), right_margin = "5"), 
          shinycssloaders::withSpinner(plotlyOutput(ns("bplot")), 
                                       type = 8, color = "green"), dataTableDownload_ui(ns("stab")))
}


exprsImpute <- function (x) 
{
  print('Running exprsImpute')
  v <- try(x@assayData$exprs_impute, silent = TRUE)
  if (inherits(v, "try-error")) 
    v <- NULL
  v
}


exprspca <- function (x, n = min(8, ncol(x) - 1), prefix = "PCA|All", fillNA = FALSE, method = 'custom',
                      ...) 
{
  print('calling exprspca function from omicsViewer')
  writePC <- function(x, n) {
    n <- min(n, length(x$sdev))
    var <- round(x$sdev[seq_len(n)]^2/(sum(x$sdev^2)) * 100, 
                 digits = 1)
    xx <- x$x[, seq_len(min(n, ncol(x$x)))]
    colnames(xx) <- paste0(prefix, "|", colnames(xx), "(", 
                           var, "%", ")")
    pp <- x$rotation[, seq_len(min(n, ncol(x$x)))]
    colnames(pp) <- paste0(prefix, "|", colnames(pp), "(", 
                           var, "%", ")")
    list(samples = xx, features = pp)
  }
  if (fillNA) {
    print(paste0(' #### using the imputation method from exprspca line 1176 omicsViewer #### ',method))
    x <- fillNA(x, method = method)
    pc <- prcomp(t(x))
  }
  else {
    nr <- nrow(x)
    x <- na.omit(x)
    pc <- prcomp(t(x), ...)
    pos <- setdiff(seq_len(nr), attr(x, "na.action"))
    rotation <- matrix(NA, nrow = nr, ncol = ncol(pc$rotation))
    rotation[pos, ] <- pc$rotation
    colnames(rotation) <- colnames(pc$rotation)
    pc$rotation <- rotation
  }
  writePC(pc, n = n)
}



#nonstandardGenericFunction for "extendMetaData" defined from package "omicsViewer"
extendMetaData <- function (object, newData, where) 
{
  standardGeneric("extendMetaData")
}


factorIndependency <- function (x, y) 
{
  tab <- table(x, y)
  r1 <- chisq.test(tab)
  r2 <- try(fisher.test(tab), silent = TRUE)
  if (is(r2, "try-error")) 
    r2 <- fisher.test(tab, simulate.p.value = TRUE, B = 1e+05)
  df <- data.frame(method = c(r1$method, r2$method), pvalue = signif(c(r1$p.value, 
                                                                       r2$p.value), digits = 3), check.names = FALSE, row.names = NULL, 
                   stringsAsFactors = FALSE)
  toMat <- function(x) apply(x, 2, function(i) i)
  if (nrow(tab) == 2 && ncol(tab) == 2) {
    odds_wald <- (tab[1, 1]/tab[1, 2])/(tab[2, 1]/tab[2, 
                                                      2])
    odds_fisher <- r2$estimate[[1]]
    df$OR <- signif(c(odds_wald, odds_fisher))
    df$OR.method <- c("Wald", "Fisher")
  }
  list(cont.table = data.frame(toMat(tab)), residual.ratio = signif(toMat(r1$observed)/toMat(r1$expected), 
                                                                    digits = 3), p.table = df)
}


factorIndependency_module <- function (input, output, session, x, y, reactive_checkpoint = reactive(TRUE)) 
{
  ns <- session$ns
  stats <- reactive({
    req(reactive_checkpoint())
    tx <- table(x())
    ty <- table(y())

    if (length(tx) < 2 || length(ty) < 2 || max(tx) < 2 || 
        max(ty) < 2 || length(tx) > 12 || length(ty) > 12) 
      return(NULL)
    factorIndependency(x = x(), y = y())
  })
  output$error <- renderUI(verbatimTextOutput(ns("errorMsg")))
  output$errorMsg <- renderText({
    req(is.null(stats()))
    "The selected variable is not suitable for independence test, too many distinct values or no duplicate value!"
  })
  output$count.table.output <- DT::renderDataTable({
    req(stats())
    DT::datatable(stats()$cont.table, options = list(searching = FALSE, 
                                                     lengthChange = FALSE, dom = "t", scrollX = TRUE), 
                  rownames = TRUE, class = "compact", caption = "Contingency table")
  })
  output$residual.ratio.output <- DT::renderDataTable({
    req(stats())
    DT::datatable(stats()$residual.ratio, options = list(searching = FALSE, 
                                                         lengthChange = FALSE, dom = "t", scrollX = TRUE), 
                  rownames = TRUE, class = "compact", caption = "Contingency table fold change (observed/expect)")
  })
  output$p.table.output <- DT::renderDataTable({
    req(stats())
    DT::datatable(stats()$p.table, options = list(searching = FALSE, 
                                                  lengthChange = FALSE, dom = "t"), rownames = FALSE, 
                  class = "compact", caption = "Significance test of the independence")
  })
}


factorIndependency_ui <- function (id) 
{
  ns <- NS(id)
  tagList(wellPanel(uiOutput(ns("error")), DT::dataTableOutput(ns("count.table.output")), 
                    DT::dataTableOutput(ns("residual.ratio.output")), DT::dataTableOutput(ns("p.table.output"))))
}


feature_general_module <- function (input, output, session, reactive_expr, reactive_i = reactive(NULL), 
                                    reactive_highlight = reactive(NULL), reactive_phenoData, 
                                    reactive_featureData, reactive_status = reactive(NULL)) 
{
  ns <- session$ns
  triset <- reactive({
    ts <- trisetter(expr = reactive_expr(), meta = reactive_phenoData(), 
                    combine = "none")
    ts[ts[, 1] != "Surv", ]
  })
  xax <- reactiveVal()
  v1 <- callModule(triselector_module, id = "tris_feature_general", 
                   reactive_x = triset, label = "Link to variable", reactive_selector1 = reactive(xax()$v1), 
                   reactive_selector2 = reactive(xax()$v2), reactive_selector3 = reactive(xax()$v3))
  attr4select_status <- reactiveVal()
  attr4select <- callModule(attr4selector_module, id = "a4_gf", 
                            reactive_meta = reactive_phenoData, reactive_expr = reactive_expr, 
                            reactive_triset = triset, reactive_status = attr4select_status)
  pheno <- reactive({
    req(v1())
    cs <- do.call(paste, list(v1(), collapse = "|"))
    if (!cs %in% colnames(reactive_phenoData())) 
      return(NULL)
    reactive_phenoData()[, cs]
  })
  pheno_cat <- reactive({
    is.factor(pheno()) || is.character(pheno())
  })
  pheno_num <- reactive({
    is.numeric(pheno())
  })
  single_i <- reactive({
    length(reactive_i()) == 1
  })
  showBoxplot <- reactive(length(pheno()) == 0 || (pheno_num() && 
                                                     !single_i()) || length(reactive_i()) == 0 || length(reactive_i()) >= 
                            10)
  showBeeswarm <- reactive(pheno_cat() && length(reactive_i()) > 
                             0 && length(reactive_i()) < 10)
  showScatter <- reactive(single_i() && pheno_num())
  output$feature_general_plot <- renderUI({
    if (showBoxplot()) 
      return(plotly_boxplot_ui(ns("feature_general_boxplotly")))
    if (showScatter()) 
      return(plotly_scatter_ui(ns("feature_general_scatter")))
    if (showBeeswarm()) {
      if (input$internal_radio == "Bees") 
        r <- plotly_scatter_ui(ns("feature_general_beeswarm"))
      else r <- plot_roc_pr_ui(ns("feature_general_roc_pr"))
      r
    }
  })
  rh <- reactiveVal()
  observeEvent(reactive_highlight(), {
    r <- reactive_highlight()
    if (is.null(r) || is.logical(r)) 
      return(NULL)
    if (length(r) == 0 && !is.null(rh()) && length(rh()) > 
        0) {
      rh(integer(0))
    }
    else if (length(r) > 0) {
      rh(match(r, colnames(reactive_expr())))
    }
  })
  callModule(plotly_boxplot_module, id = "feature_general_boxplotly", 
             reactive_param_plotly_boxplot = reactive({
               req(reactive_expr())
               ylab <- rownames(reactive_expr())[reactive_i()]
               if (length(ylab) > 1) 
                 ylab <- "Abundance of selected features"
               else if (length(ylab) == 0) 
                 ylab <- "Relative abundance"
               else if (is.na(ylab)) 
                 ylab <- "Relative abundance"
               ylab.extvar <- do.call(paste, list(v1(), collapse = "|"))
               list(x = reactive_expr(), i = reactive_i(), highlight = rh(), 
                    extvar = pheno(), ylab = ylab, ylab.extvar = ylab.extvar)
             }), reactive_checkpoint = showBoxplot)
  scatter_vars <- reactive({
    l <- list(source = "feature_general_module")
    l$color <- attr4select$color
    l$shape <- attr4select$shape
    l$size <- attr4select$size
    l$tooltips <- attr4select$tooltips
    l$highlight <- attr4select$highlight
    l$highlightName <- attr4select$highlightName
    req(reactive_i() %in% rownames(reactive_expr()) || reactive_i() <= 
          nrow(reactive_expr()))
    if (showScatter()) {
      l$x <- reactive_expr()[reactive_i(), ]
      l$y <- pheno()
      l$xlab <- rownames(reactive_expr())[reactive_i()]
      l$ylab <- do.call(paste, list(v1(), collapse = "|"))
      if (is.null(l$tooltips)) 
        l$tooltips <- colnames(reactive_expr())
    }
    if (showBeeswarm()) {
      df <- reshape2::melt(reactive_expr()[reactive_i(), , drop = FALSE])
      df$color <- rep(l$color, each = length(reactive_i()))
      df$pheno <- rep(pheno(), each = length(reactive_i()))
      xlab <- ""
      ylab <- rownames(reactive_expr())[reactive_i()]
      if (length(ylab) > 1) 
        ylab <- "Abundance of selected features"
      else if (length(ylab) == 0) 
        ylab <- "Relative abundance"
      else if (is.na(ylab)) 
        ylab <- "Relative abundance"
      df <- na.omit(df)
      l$y <- df$value
      l$x <- df$pheno
      l$ylab <- ylab
      l$color <- df$color
      if (is.null(l$tooltips)) 
        l$tooltips <- sprintf("<b>Feature: </b>%s<br><b>Sample: </b>%s", 
                              df$Var1, df$Var2)
    }
    l
  })
  showRegLine <- reactiveVal(FALSE)
  htestV1 <- reactiveVal()
  htestV2 <- reactiveVal()
  v_scatter <- callModule(plotly_scatter_module, id = "feature_general_scatter", 
                          reactive_param_plotly_scatter = scatter_vars, reactive_checkpoint = showScatter, 
                          reactive_regLine = reactive(showRegLine()))
  observe({
    showRegLine(v_scatter()$regline)
  })
  v_beeswarm <- callModule(plotly_scatter_module, id = "feature_general_beeswarm", 
                           reactive_param_plotly_scatter = scatter_vars, reactive_checkpoint = showBeeswarm, 
                           htest_var1 = htestV1, htest_var2 = htestV2)
  callModule(plot_roc_pr_module, id = "feature_general_roc_pr", 
             reactive_param = scatter_vars, reactive_checkpoint = showBeeswarm)
  metatab <- reactive({
    req(reactive_i())
    tab <- reactive_featureData()
    tab <- tab[, grep("^General\\|", colnames(tab)), drop = FALSE]
    tab <- tab[reactive_i(), , drop = FALSE]
    ic <- vapply(tab, is.numeric, logical(1)) & vapply(tab, 
                                                       is.integer, logical(1))
    # tab[ic] <- lapply(tab[ic], signif, digits = 2)
    colnames(tab) <- sub("General\\|All\\|", "", colnames(tab))
    tab
  })
  callModule(dataTableDownload_module, id = "mtab", reactive_table = metatab, 
             prefix = "FeatureTable_")
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status())) 
      return()
    xax(NULL)
    xax(list(v1 = s$xax[[1]], v2 = s$xax[[2]], v3 = s$xax[[3]]))
  })
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status())) 
      return()
    attr4select_status(NULL)
    attr4select_status(s$attr4)
  })
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status())) 
      return()
    htestV1(s$htestV1)
    htestV2(s$htestV2)
  })
  observeEvent(reactive_status(), {
    if (!is.null(s <- reactive_status())) 
      showRegLine(s$showRegLine)
  })
  rv <- reactiveValues()
  observe(rv$xax <- v1())
  observe(rv$showRegLine <- showRegLine())
  observe(rv$attr4 <- attr4select$status)
  observe({
    rv$htestV1 <- v_beeswarm()$htestV1
    rv$htestV2 <- v_beeswarm()$htestV2
  })
  reactive({
    reactiveValuesToList(rv)
  })
}


feature_general_ui <- function (id) 
{
  ns <- NS(id)
  tagList(fluidRow(column(12, style = "margin-top: 0px;", triselector_ui(ns("tris_feature_general"), 
                                                                         right_margin = "5")), column(11, uiOutput(ns("feature_general_plot"))), 
                   column(1, attr4selector_ui(ns("a4_gf"), circle = FALSE, 
                                              right = TRUE), radioGroupButtons(inputId = ns("internal_radio"), 
                                                                               label = " ", size = "xs", choices = c("Bees", "Curve"), 
                                                                               direction = "vertical", status = "success"))), dataTableDownload_ui(ns("mtab")))
}


fgsea1 <- function (gs, stats, gs_desc = NULL, ...) 
{
  if (inherits(gs, c("matrix", "dgCMatrix"))) {
    gs <- csc2list(gs)
  }
  pw <- split(gs$featureId, gs$gsId)
  params <- list(pathways = pw, stats = stats, ...)
  res <- do.call(fgseaMultilevel, params)
  if (!is.null(gs_desc)) {
    res$desc <- as.character(gs_desc[res$pathway])
  }
  as.data.frame(res)
}



fillNA <- function(x, method='perseus'){
  result = x
  result <- tryCatch({
  if ((method == 'custom') | (method == 'chen-meng')) {
    result = impute_custom(x)
  } else {
    result = impute_perseus(x)
  }
}, error = function(e) {
  message("Error caught: ", e$message)
  result  # Return fallback value
})
  print(any(is.na(result)))
  print(mean(result))
  return(result)
} 


impute_custom <- function (x, maxfill = quantile(x, probs = 0.15, na.rm = TRUE), 
                    fillingFun = function(x) min(x, na.rm = TRUE) - log10(2)) 
{
  print('impute_custom: Running imputation function for Chen method')
  xf <- apply(x, 1, function(xx) {
    x3 <- xx
    x3[is.na(x3)] <- min(maxfill, fillingFun(xx))
    x3
  })
  xf <- t(xf)
  rownames(xf) <- rownames(x)
  colnames(xf) <- colnames(xf)
  xf
}


impute_perseus <- function(object, width=0.3, downshift=1.8, seed=100) {
  print('impute_perseus: Running imputation function for perseus method')
  mx <- max(object, na.rm=TRUE)
  mn <- min(object, na.rm=TRUE)
  set.seed(seed)
  xf <- apply(object, 2, function(temp) {
    temp[!is.finite(temp)] <- NA
    temp_sd <- stats::sd(temp, na.rm=TRUE)
    temp_mean <- mean(temp, na.rm=TRUE)
    shrinked_sd <- width * temp_sd   # shrink sd width
    downshifted_mean <- temp_mean - downshift * temp_sd   # shift mean of imputed values
    n_missing <- sum(is.na(temp))
    temp[is.na(temp)] <- stats::rnorm(n_missing, mean=downshifted_mean, sd=shrinked_sd)
    temp
  })

  xf
}





filterRow <- function (x, max.quantile = NULL, max.value = NULL, var = NULL, 
                       min.rep = 2) 
{
  rmi <- rowMaxs(x, na.rm = TRUE)
  if (!is.null(max.value)) {
    qt <- max.value
  }
  else if (!is.null(max.quantile)) {
    qt <- quantile(x, na.rm = TRUE, probs = max.quantile)
  }
  else qt <- NULL
  f1 <- FALSE
  if (!is.null(qt)) 
    f1 <- rmi < qt
  f2 <- FALSE
  if (!is.null(var)) {
    sn <- vapply(unique(var), function(xx) {
      rowSums(!is.na(x[, var == xx, drop = FALSE]))
    }, double(nrow(x)))
    rsn <- rowMaxs(sn)
    f2 <- rsn < min.rep
  }
  !(f1 | f2)
}


geneshot_module <- function (input, output, session, pdata, fdata, expr, feature_selected, 
                             sample_selected, object, reactive_status = reactive(NULL)) 
{
  ns <- session$ns
  triset <- reactive({
    fd <- fdata()
    cn <- colnames(fd)[!vapply(fd, is.numeric, logical(1)) & 
                         !grepl("^GS\\|", colnames(fd))]
    str_split_fixed(cn, "\\|", n = 3)
  })
  xax <- reactiveVal()
  v1 <- callModule(triselector_module, id = "geneNameCol", 
                   reactive_x = triset, label = "Map ID", reactive_selector1 = reactive(xax()$v1), 
                   reactive_selector2 = reactive(xax()$v2), reactive_selector3 = reactive(xax()$v3))
  rif <- reactiveVal()
  observeEvent(input$submit, {
    show_modal_spinner(text = "Querying database ...")
    res <- getAutoRIF(trimws(strsplit(input$term, ";")[[1]]), 
                      filter = TRUE)
    if (!is.null(res) && nrow(res) > 0) {
      res$selected <- ""
      res <- res[order(res$rank, decreasing = TRUE), ]
      dft <- res[seq_len(min(nrow(res), 20)), ]
      outliers <- list(x = dft$n, y = dft$perc, text = dft$gene, 
                       xref = "x", yref = "y", showarrow = TRUE, ax = 10, 
                       ay = -20)
    }
    remove_modal_spinner()
    if (is.null(res) || nrow(res) == 0) 
      showModal(modalDialog("No genes identified, perhaps the terms is not correctly typed."))
    req(!is.null(res) && nrow(res) > 0)
    rif(list(tab = res, outliers = outliers))
  })
  outliersLabs <- reactive(rif()$outliers)
  gtab <- reactive({
    req(df <- rif()$tab)
    cid <- paste(unlist(v1()), collapse = "|")
    if (cid %in% colnames(fdata())) 
      df$selected[df$gene %in% fdata()[feature_selected(), 
                                       cid]] <- "+"
    df
  })
  rifRow <- callModule(dataTableDownload_module, id = "autorif", 
                       reactive_table = reactive({
                         req(gtab())
                         gtab()
                       }), prefix = "autoRIF", pageLength = 10)
  output$plt <- plotly::renderPlotly({
    req(df <- gtab())
    fig <- plotly_scatter(x = df$n, y = df$perc, xlab = "# of publication", 
                          ylab = "Publications with Search Term(s) / Total Publications", 
                          color = df$selected, size = 10, tooltips = df$gene, 
                          shape = "select")
    fig <- plotly::layout(fig$fig, annotations = outliersLabs())
    plotly::config(fig, toImageButtonOptions = list(format = "svg", 
                                                    filename = "omicsViewerPlot", width = 700, height = 700))
  })
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status())) 
      return()
    xax(NULL)
    xax(list(v1 = s$xax[[1]], v2 = s$xax[[2]], v3 = s$xax[[3]]))
    rif(s$rif)
    updateTextInput(session, inputId = "term", value = s$term)
  })
  rv <- reactiveValues()
  observe(rv$xax <- v1())
  observe(rv$rif <- rif())
  observe(rv$term <- input$term)
  reactive({
    reactiveValuesToList(rv)
  })
}


geneshot_ui <- function (id) 
{
  ns <- NS(id)
  tagList(fluidRow(column(width = 8, textInput(ns("term"), 
                                               label = NULL, value = "", width = "100%", placeholder = "Search any term here, multiple separated by (;) ...")), 
                   column(width = 4, align = "right", actionButton(ns("submit"), 
                                                                   label = "Search related genes!", width = "100%")), 
                   column(width = 11, triselector_ui(ns("geneNameCol")))), 
          plotly::plotlyOutput(ns("plt")), tags$br(), dataTableDownload_ui(ns("autorif")))
}


getAutoRIF <- function (term, rif = c("generif", "autorif")[1], filter = TRUE) 
{
  term <- gsub(" ", "%20", term)
  term <- paste(term, collapse = ",")
  GENESHOT_URL <- "https://maayanlab.cloud/geneshot/api/search"
  payload <- list(rif = rif, term = term)
  r <- httr::POST(GENESHOT_URL, body = payload, encode = "json")
  r <- httr::content(r)
  v <- vapply(r$gene_count, unlist, numeric(2))
  if (length(v) == 0) 
    return(NULL)
  df <- data.frame(gene = colnames(v), n = v[1, ], perc = v[2, 
  ], stringsAsFactors = FALSE)
  df$rank <- df$n * df$perc
  if (filter) 
    df <- df[df$n > min(ceiling(nrow(df)/200), 3), ]
  attr(df, "term") <- r$search_term
  attr(df, "pubmedID_count") <- r$pubmedID_count
  df
}


getAx <- function (x, what) 
{
  if (inherits(x, "SQLiteConnection")) {
    if (what == "dendrogram") {
      v <- getDend(x)
    }
    else {
      v <- dbGetQuery(x, "SELECT value FROM axes WHERE axis = :x;", 
                      params = list(x = what))[[1]]
      if (length(v) == 0) 
        v <- NULL
    }
  }
  else if (inherits(x, "ExpressionSet")) {
    v <- attr(x, what)
  }
  v
}


getDend <- function (x) 
{
  dend <- dbGetQuery(x, "SELECT * FROM dendrogram;")
  if (nrow(dend) == 0) 
    return(NULL)
  l <- lapply(dend$object, function(x) {
    o <- str2hclust(x)
    list(ord = o$order, hcl = as.dendrogram(o))
  })
  names(l) <- dend$name
  l
}


getExprs <- function (x) 
{
  if (inherits(x, "SQLiteConnection")) {
    mat <- dbGetQuery(x, "SELECT * FROM exprs;")
    rn <- mat$rowname
    mat$rowname <- NULL
    mat <- apply(mat, 2, as.numeric)
    rownames(mat) <- rn
  }
  else if (inherits(x, "ExpressionSet")) 
    mat <- Biobase::exprs(x)
  mat
}


getExprsImpute <- function (x) 
{
  print('Running getExprsImpute')
  if (inherits(x, "SQLiteConnection")) {
    if (!"exprsimpute" %in% dbListTables(x)) 
      return(NULL)
    mat <- dbGetQuery(x, "SELECT * FROM exprsimpute;")
    rn <- mat$rowname
    mat$rowname <- NULL
    mat <- apply(mat, 2, as.numeric)
    rownames(mat) <- rn
  }
  else if (inherits(x, "ExpressionSet")) 
    mat <- exprsImpute(x)
  else mat <- NULL
  mat
}


getFData <- function (x) 
{
  if (inherits(x, "SQLiteConnection")) {
    mat <- dbGetQuery(x, "SELECT * FROM feature;")
    rownames(mat) <- mat$rowname
    mat$rowname <- NULL
    gs <- dbGetQuery(x, "SELECT * FROM GS;")
    if (nrow(gs) > 0) {
      gs$featureId <- as.factor(gs$featureId)
      gs$gsId <- as.factor(gs$gsId)
      attr(mat, "GS") <- gs
    }
  }
  else if (inherits(x, "ExpressionSet")) {
    mat <- fData(x)
  }
  mat
}

getMQParams <- function (x) 
{
  mqpar <- fxml_importXMLFlat(x)
  method <- "LF"
  if (length(grep("TMT[0-9]*plex", mqpar$value.)) > 1) 
    method <- "TMT"
  ifa <- which(mqpar$elem. == "fastaFilePath")
  if (length(ifa) == 0) 
    ifa <- which(mqpar$elem. == "string" & mqpar$level2 == 
                   "fastaFiles")
  list(version = mqpar$value.[which(mqpar$elem. == "maxQuantVersion")], 
       fasta = mqpar$value.[ifa], enzymes = na.omit(mqpar$value.[which(mqpar$level4 == 
                                                                         "enzymes")]), varMod = na.omit(mqpar$value.[which(mqpar$level4 == 
                                                                                                                             "variableModifications")]), fixedMod = na.omit(mqpar$value.[which(mqpar$level4 == 
                                                                                                                                                                                                 "fixedModifications")]), mainTol = mqpar$value.[which(mqpar$elem. == 
                                                                                                                                                                                                                                                         "mainSearchTol")], MBR = mqpar$value.[which(mqpar$elem. == 
                                                                                                                                                                                                                                                                                                       "matchBetweenRuns")], peptideFdr = mqpar$value.[which(mqpar$elem. == 
                                                                                                                                                                                                                                                                                                                                                               "peptideFdr")], proteinFdr = mqpar$value.[which(mqpar$elem. == 
                                                                                                                                                                                                                                                                                                                                                                                                                 "proteinFdr")], ibaq = mqpar$value.[which(mqpar$elem. == 
                                                                                                                                                                                                                                                                                                                                                                                                                                                             "ibaq")], lfqMode = mqpar$value.[which(mqpar$elem. == 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "lfqMode")], method = method)
}


getOrderCols <- function (x) 
{
  if (is.null(x)) 
    return(NULL)
  x[["order"]]
}


getOutliersDown <- function (x, threshold = 0.1) 
{
  f <- apply(x, 1, function(x) {
    x <- sort(x, decreasing = FALSE)
    if (is.na(x[2])) 
      return(c(NA, NA))
    c(min(x[1] - x[2], 0), names(x[1]))
  })
  df <- data.frame(fold.change.log10 = as.numeric(f[1, ]), 
                   sample = f[2, ], stringsAsFactors = FALSE)
  df[which(abs(df$fold.change.log10) < threshold), ] <- NA
  df
}

getOutliersUp <- function (x, na.replace = quantile(x, 0.25, na.rm = TRUE), threshold = 0.1) 
{
  sec <- na.replace
  f <- apply(x, 1, function(x, sec) {
    x <- sort(x, decreasing = TRUE)
    if (!is.na(x[2])) 
      sec <- x[2]
    c(max(x[1] - sec, 0), names(x[1]))
  }, sec = sec)
  df <- data.frame(fold.change.log10 = as.numeric(f[1, ]), 
                   sample = f[2, ], stringsAsFactors = FALSE)
  df[which(abs(df$fold.change.log10) < threshold), ] <- NA
  df
}


getPData <- function (x) 
{
  if (inherits(x, "SQLiteConnection")) {
    mat <- dbGetQuery(x, "SELECT * FROM sample;")
    rownames(mat) <- mat$rowname
    mat$rowname <- NULL
  }
  else if (inherits(x, "ExpressionSet")) {
    mat <- pData(x)
  }
  mat
}


getSearchCols <- function (x) 
{
  if (is.null(x)) 
    return(NULL)
  lapply(x$columns, "[[", "search")
}


getSteps <- function (fdata, stepList) 
{
  sts <- lapply(stepList, function(cc) {
    cc <- unlist(cc)
    pname <- cc[1]
    gname <- cc[-1]
    n <- length(gname)
    m1 <- fdata[[paste("mean", pname, gname[1], sep = "|")]]
    mn <- fdata[[paste("mean", pname, gname[n], sep = "|")]]
    mm <- abs(m1 - mn)
    fillFalse <- function(x) {
      x[is.na(x)] <- FALSE
      x
    }
    s <- lapply(seq_len(length(gname) - 1), function(i) {
      x <- gname[c(i, i + 1)]
      ifdr <- sprintf("ttest|%s_vs_%s|fdr", x[1], x[2])
      imd <- sprintf("ttest|%s_vs_%s|mean.diff", x[1], 
                     x[2])
      list(signif = cbind(fdata[, ifdr] < 0.01, fdata[, 
                                                      ifdr] < 0.05, fdata[, ifdr] < 0.1), down = fillFalse(fdata[, 
                                                                                                                 imd] > 0), up = fillFalse(fdata[, imd] < 0))
    })
    down <- lapply(s, "[[", "down")
    down <- Reduce("&", down)
    up <- lapply(s, "[[", "up")
    up <- Reduce("&", up)
    sig <- lapply(s, "[[", "signif")
    sig <- Reduce("+", sig)
    sig[!(up | down), ] <- 0
    sig[down, ] <- -sig[down, ]
    colnames(sig) <- paste("x.fdr", c("01", "05", "1"), sep = ".")
    res <- cbind(sig, y.mean.diff = mm)
    colnames(res) <- paste(paste(cc, collapse = ":"), colnames(res), 
                           sep = "|")
    res
  })
  r <- do.call(cbind, sts)
  colnames(r) <- paste("step", colnames(r), sep = "|")
  r
}


getStringId <- function (genes, taxid = 9606, caller = "omicsViewer") 
{
  string_api_url <- "https://string-db.org/api"
  output_format <- "tsv"
  method <- "get_string_ids"
  params <- list(identifiers = paste(genes, collapse = "%0d"), 
                 species = taxid, limit = 1, echo_query = 1, caller_identity = caller)
  params <- paste(vapply(names(params), function(x) paste(x, 
                                                          params[[x]], sep = "="), FUN.VALUE = character(1)), collapse = "&")
  request_url <- paste(string_api_url, output_format, method, 
                       sep = "/")
  results <- httr::GET(paste(request_url, params, sep = "?"))
  httr::content(results)
}


getUPRefProteomeID <- function (domain = c("Eukaryota", "Archaea", "Bacteria", "Viruses")[1]) 
{
  url <- curl(sprintf("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/%s/", 
                      domain))
  on.exit(close(url))
  con <- readLines(url)
  na.omit(stringr::str_match(con, "\\>\\s*(.*?)\\s*/\\<")[, 
                                                          2])
}


gsAnnotIdList <- function (idList, gsIdMap, minSize = 5, maxSize = 500, data.frame = FALSE, 
                           sparse = TRUE) 
{
  pid <- data.frame(id = unlist(idList), index = rep(seq_along(idList), 
                                                     times = vapply(idList, length, integer(1))), stringsAsFactors = FALSE)
  gsIdMap <- gsIdMap[gsIdMap$id %in% pid$id, ]
  if (data.frame) {
    v <- split(pid$index, pid$id)[gsIdMap$id]
    v <- data.frame(featureId = seq_along(idList)[unlist(v)], 
                    gsId = rep(gsIdMap$term, vapply(v, length, integer(1))), 
                    weight = 1, stringsAsFactors = FALSE)
    v <- unique(v)
    n <- table(v$gsId)
    exc <- names(n[n > maxSize | n < minSize])
    v <- v[!v$gsId %fin% exc, ]
  }
  else {
    gsIdMap <- split(gsIdMap$id, gsIdMap$term)
    v <- vapply(gsIdMap, function(x, s) {
      s[unique(pid$index[fastmatch::"%fin%"(pid$id, x)])] <- 1
      s
    }, s = rep(0, length(idList)), FUN.VALUE = numeric(length(idList)))
    colnames(v) <- names(gsIdMap)
    nm <- matrixStats::colSums2(v)
    v <- v[, which(nm >= minSize & nm <= maxSize)]
    if (sparse) 
      v <- as(v, "dgCMatrix")
  }
  v
}


gslist_module <- function (input, output, session, reactive_featureData, reactive_i) 
{
  ns <- session$ns
  reactive_pathway <- reactive({
    req(f1 <- reactive_featureData())
    gss <- attr(f1, "GS")
    req(gss)
    s <- cbind(gss, f1[fmatch(gss$featureId, rownames(f1)), 
                       grep("^General", colnames(f1)), drop = FALSE])
    colnames(s)[colnames(s) == "gsId"] <- "Gene-set"
    s
  })
  tab <- reactive({
    req(reactive_pathway())
    if (length(reactive_i()) == 1 && is.logical(reactive_i()) && reactive_i())   return(reactive_pathway())
    if (length(reactive_i()) == 0 || all(is.na(reactive_i())))   return(reactive_pathway())

    df <- reactive_pathway()[reactive_pathway()$featureId %fin%  reactive_i(), ]
    req(is.data.frame(df))
    df
  })
  ii <- callModule(dataTableDownload_module, id = "stab", reactive_table = reactive({
    tab()[, setdiff(colnames(tab()), "featureId")]
  }), prefix = "gslist_", pageLength = 25)
  reactive({
    req(ii())
    as.character(tab()$featureId[ii()])
  })
}


gslist_ui <- function (id) 
{
  ns <- NS(id)
  dataTableDownload_ui(ns("stab"))
}


hasAttr <- function (x, attr.name) 
{
  attr.name %in% names(attributes(x))
}


hclust2str <- function (x) 
{
  cc <- c(merge = paste(paste(x$merge[, 1], x$merge[, 2], sep = "\t"), 
                        collapse = "\n"), height = paste(signif(x$height, digits = 5), 
                                                         collapse = ";"), order = paste(x$order, collapse = ";"), 
          labels = paste(x$labels, collapse = "=;="), method = x$method, 
          dist.method = x$dist.method)
  paste(cc, collapse = "__elementSplitter__")
}


heatmapKey <- function (range, colors) 
{
  par(mar = c(2, 1, 0, 1))
  bp <- barplot(rep(1, length(colors)), col = colors, space = 0, 
                axes = FALSE, xaxs = "i", yaxs = "i", border = colors)
  lab <- grid.pretty(range)
  at <- (lab - min(range))/(max(range) - min(range)) * max(bp)
  mtext(text = lab, side = 1, at = at)
}


# iheatmap <- function (x, fData = NULL, pData = NULL, impute = FALSE) 
# {
#   if (inherits(x, "ExpressionSet") || inherits(x, "xcmsFeatureSet")) {
#     fData <- fData(x)
#     pData <- pData(x)
#     x <- Biobase::exprs(x)
#   }
#   ir <- unique(c(which(rowSums2(!is.na(x)) == 0), which(rowVars(x) == 
#                                                           0)))
#   if (length(ir) > 0) {
#     fData <- fData[-ir, ]
#     x <- x[-ir, ]
#   }
#   if (impute) {
#     x <- apply(x, 1, function(xx) {
#       xx[is.na(xx)] <- min(xx, na.rm = TRUE) * 0.9
#       xx
#     })
#     x <- t(x)
#   }
#   ui <- fluidPage(sidebarLayout(sidebarPanel = sidebarPanel(tabsetPanel(tabPanel("Parameters", 
#                               iheatmapInput(id = "test")), tabPanel("Legend", iheatmapLegend(id = "test")))), 
#                                 mainPanel = mainPanel(iheatmapOutput(id = "test"))))

#   server <- function(input, output) {

#     callModule(iheatmapModule,
#                "test",
#                mat = reactive(x), 
#                pd = reactive(pData),
#                fd = reactive(fData))

#   }
#   shinyApp(ui, server)
# }


iheatmapClear <- function (id) 
{
  ns <- NS(id)
  actionBttn(ns("clear"), "Clear figure selection", style = "minimal", 
             color = "primary", size = "xs")
}


iheatmapInput <- function (id, scaleOn = "row") 
{
  ns <- NS(id)
  tagList(fluidRow(column(6, selectInput(ns("colSortBy"), "Sorting columns by", 
                                         choices = NULL, selectize = TRUE), conditionalPanel(sprintf("input['%s'] == 'hierarchical cluster'", 
                                                                                                     ns("colSortBy")), selectInput(ns("clusterColDist"), "Distance", 
                                                                                                                                   choices = c("Pearson correlation", "Euclidean", "Maximum", 
                                                                                                                                               "Manhattan", "Canberra", "Binary", "Minkowski", "Spearman correlation"), 
                                                                                                                                   selectize = TRUE), selectInput(ns("clusterColLink"), 
                                                                                                                                                                  "Linkage", choices = c("ward.D", "ward.D2", "single", 
                                                                                                                                                                                         "complete", "average", "mcquitty", "median", "centroid"), 
                                                                                                                                                                  selectize = TRUE)), hr(), selectInput(ns("rowSortBy"), 
                                                                                                                                                                                                        "Sorting rows by", choices = NULL, selectize = TRUE, 
                                                                                                                                                                                                        multiple = FALSE), conditionalPanel(sprintf("input['%s'] == 'hierarchical cluster'", 
                                                                                                                                                                                                                                                    ns("rowSortBy")), selectInput(ns("clusterRowDist"), "Distance", 
                                                                                                                                                                                                                                                                                  choices = c("Pearson correlation", "Euclidean", "Maximum", 
                                                                                                                                                                                                                                                                                              "Manhattan", "Canberra", "Binary", "Minkowski", "Spearman correlation"), 
                                                                                                                                                                                                                                                                                  selectize = TRUE), selectInput(ns("clusterRowLink"), 
                                                                                                                                                                                                                                                                                                                 "Linkage", choices = c("ward.D", "ward.D2", "single", 
                                                                                                                                                                                                                                                                                                                                        "complete", "average", "mcquitty", "median", "centroid"), 
                                                                                                                                                                                                                                                                                                                 selectize = TRUE)), hr(), selectInput(ns("scale"), "Scale on", 
                                                                                                                                                                                                                                                                                                                                                       choices = c("row", "none", "column"), selected = scaleOn, 
                                                                                                                                                                                                                                                                                                                                                       selectize = TRUE)), column(6, selectInput(ns("annotCol"), 
                                                                                                                                                                                                                                                                                                                                                                                                 label = "Column annotations", choices = NULL, multiple = TRUE), 
                                                                                                                                                                                                                                                                                                                                                                                  selectInput(ns("annotRow"), label = "Row annotations", 
                                                                                                                                                                                                                                                                                                                                                                                              choices = NULL, multiple = TRUE), selectInput(ns("tooltipInfo"), 
                                                                                                                                                                                                                                                                                                                                                                                                                                            label = "Tooltips", choices = NULL, multiple = TRUE), 
                                                                                                                                                                                                                                                                                                                                                                                  hr(), sliderInput(ns("marginBottom"), "Bottom margin", 
                                                                                                                                                                                                                                                                                                                                                                                                    min = 1, max = 20, value = 4), sliderInput(ns("marginRight"), 
                                                                                                                                                                                                                                                                                                                                                                                                                                               "Right margin", min = 1, max = 20, value = 4), selectInput(ns("heatmapColors"), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          label = "Heatmap color panel", choices = c("BrBG", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "RdYlGn"), selected = "RdYlBu", selectize = TRUE, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          multiple = FALSE))))
}


iheatmapLegend <- function (id) 
{
  ns <- NS(id)
  tagList(uiOutput(ns("key_heatmap_ui")), uiOutput(ns("key_colSideCor_ui")), 
          uiOutput(ns("key_rowSideCor_ui")))
}


iheatmapModule <- function (input, output, session, mat, pd, fd, status = reactive(NULL), 
                            fill.NA = TRUE , method = 'perseus') 
{
  print('Running iheatmapModule')
  ns <- session$ns
  matr <- reactive({
    req(mat())
    if (any(is.na(mat())) && fill.NA) 
      r <- fillNA(mat(),method = method)
    else r <- mat()
    r
  })
  rowDendrogram <- reactive({
    req(mat())
    attr(mat(), "rowDendrogram")
  })
  colDendrogram <- reactive({
    req(mat())
    attr(mat(), "colDendrogram")
  })
  clsRow <- reactive({
    if (nrow(matr()) > 750 || !is.null(rowDendrogram())) 
      return("none")
    "hierarchical cluster"
  })
  rdg <- reactive({
    if (is.null(rowDendrogram()) || is.null(names(rowDendrogram()))) 
      return(NULL)
    x <- rowDendrogram()
    names(x) <- paste("HCL", names(x))
    x
  })
  cdg <- reactive({
    if (is.null(colDendrogram()) || is.null(names(colDendrogram()))) 
      return(NULL)
    x <- colDendrogram()
    names(x) <- paste("HCL", names(x))
    x
  })
  gsdf <- reactive({
    attr(fd(), "GS")
  })
  fdColWithGS <- reactive({
    req(fd())
    c(colnames(fd()), paste0("GS|", levels(gsdf()$gsId)))
  })
  observe({
    req(pd())
    updateSelectInput(session, "annotCol", choices = colnames(pd()))
  })
  observe({
    req(pd())
    updateSelectInput(session, "annotRow", choices = fdColWithGS())
  })
  observe({
    req(pd())
    updateSelectInput(session, "colSortBy", choices = c(names(cdg()), 
                                                        "hierarchical cluster", "none", colnames(pd())))
  })
  observe({
    cs <- c(names(rdg()), "none", "hierarchical cluster", 
            fdColWithGS())
    ss <- clsRow()
    if (ss == "none" && cs[[1]] != "none") 
      ss <- cs[[1]]
    updateSelectInput(session, "rowSortBy", choices = cs, 
                      selected = ss)
  })
  observe({
    req(pd())
    updateSelectInput(session, "tooltipInfo", choices = c(colnames(fd()), 
                                                          colnames(pd())))
  })
  mm <- reactive({
    req(input$scale)
    req(matr())
    if (input$scale == "row") {
      mm <- t(scale(t(matr())))
      brk <- c(min(mm, na.rm = TRUE), seq(-2, 2, length.out = 99), 
               max(mm, na.rm = TRUE))
    }
    else if (input$scale == "column") {
      mm <- scale(matr())
      brk <- c(min(mm, na.rm = TRUE), seq(-2, 2, length.out = 99), 
               max(mm, na.rm = TRUE))
    }
    else {
      mm <- matr()
      brk <- seq(min(mm, na.rm = TRUE), max(mm, na.rm = TRUE), 
                 length.out = 101)
    }
    list(mat = mm, breaks = brk)
  })
  pre_hcl <- reactiveVal()
  pre_ord <- reactiveVal()
  pre_hcl_col <- reactiveVal()
  pre_ord_col <- reactiveVal()
  observeEvent(status(), {
    if (is.null(status())) 
      return(NULL)
    updateSelectInput(session, "annotCol", selected = null2empty(status()$annotCol))
    updateSelectInput(session, "annotRow", selected = null2empty(status()$annotRow))
    updateSelectInput(session, "colSortBy", selected = status()$colSortBy)
    updateSelectInput(session, "rowSortBy", selected = status()$rowSortBy)
    updateSelectInput(session, "tooltipInfo", selected = null2empty(status()$tooltipInfo))
    updateSelectInput(session, "heatmapColors", selected = status()$heatmapColors)
    updateSelectInput(session, "scale", selected = status()$scale)
    updateSelectInput(session, "clusterColDist", selected = status()$clusterColDist)
    updateSelectInput(session, "clusterColLink", selected = status()$clusterColLink)
    updateSelectInput(session, "clusterRowDist", selected = status()$clusterRowDist)
    updateSelectInput(session, "clusterRowLink", selected = status()$clusterRowLink)
    updateSliderInput(session, "marginRight", value = status()$marginRight)
    updateSliderInput(session, "marginBottom", value = status()$marginBottom)
    pre_hcl(status()$rowDendrogram)
    pre_ord(status()$rowOrder)
    pre_hcl_col(status()$colDendrogram)
    pre_ord_col(status()$colOrder)
  })
  rowSB <- eventReactive(list(input$rowSortBy, mm()$mat, input$clusterRowDist, 
                              input$clusterRowLink, status()), {
                                if (!is.null(pre_hcl()) & !is.null(pre_ord())) {
                                  return(list(ord = pre_ord(), hcl = pre_hcl()))
                                }
                                req(input$rowSortBy)
                                hcl_r <- NULL
                                ord_r <- seq_len(nrow(mm()$mat))
                                if (input$rowSortBy %in% names(rdg())) {
                                  return(rdg()[[input$rowSortBy]])
                                }
                                else if (!is.null(gsdf()) && grepl("GS\\|", input$rowSortBy)[1]) {
                                  d <- gsdf()$featureId[gsdf()$gsId == sub("^GS\\|", 
                                                                           "", input$rowSortBy)]
                                  d <- as.numeric(rownames(fd()) %fin% as.character(d))
                                  ord_r <- order(d)
                                }
                                else if (!input$rowSortBy %in% c("", "none", "hierarchical cluster")) {
                                  req(input$rowSortBy %in% colnames(fd()))
                                  ord_r <- order(fd()[, input$rowSortBy])
                                }
                                else if (input$rowSortBy == "hierarchical cluster" && 
                                         nrow(mm()$mat) > 2) {
                                  dd <- tolower(strsplit(input$clusterRowDist, " ")[[1]][1])
                                  hcl_r <- hclust(adist(mm()$mat, method = dd), method = input$clusterRowLink)
                                  ord_r <- hcl_r$order
                                  hcl_r <- as.dendrogram(hcl_r)
                                }
                                list(ord = ord_r, hcl = hcl_r)
                              })
  colSB <- eventReactive(list(input$colSortBy, mm()$mat, input$clusterColDist, 
                              input$clusterColLink, status()), {
                                if (!is.null(pre_hcl_col()) & !is.null(pre_ord_col())) {
                                  return(list(ord = pre_ord_col(), hcl = pre_hcl_col()))
                                }
                                req(input$colSortBy)
                                hcl_c <- NULL
                                ord_c <- seq_len(ncol(mm()$mat))
                                if (input$colSortBy %in% names(cdg())) {
                                  return(cdg()[[input$colSortBy]])
                                }
                                else if (!input$colSortBy %in% c("", "none", "hierarchical cluster")) {
                                  req(input$colSortBy %in% colnames(pd()))
                                  ord_c <- order(pd()[, input$colSortBy])
                                }
                                else if (input$colSortBy == "hierarchical cluster" && 
                                         ncol(mm()$mat) > 2) {
                                  dd <- tolower(strsplit(input$clusterColDist, " ")[[1]][1])
                                  hcl_c <- hclust(adist(t(mm()$mat), method = dd), 
                                                  method = input$clusterColLink)
                                  ord_c <- hcl_c$order
                                  hcl_c <- as.dendrogram(hcl_c)
                                }
                                list(ord = ord_c, hcl = hcl_c)
                              })
  hm <- reactive({
    mmo <- mm()$mat[rowSB()$ord, colSB()$ord]
    list(dend_c = colSB()$hcl, dend_r = rowSB()$hcl, mat = t(mmo), 
         ord_c = colSB()$ord, ord_r = rowSB()$ord, brk = mm()$breaks)
  })
  dat_colSideCol <- reactive({
    req(input$annotCol)
    addHeatmapAnnotation(pd()[hm()$ord_c, input$annotCol], 
                         var.name = input$annotCol)
  })
  output$colSideCol <- renderPlot({
    par(mar = c(0, 0, 0, input$marginRight))
    addHeatmapAnnotation_plot(dat_colSideCol(), xlim = ranges$x - 
                                0.5)
  })
  dat_rowSideCol <- reactive({
    req(input$annotRow)
    ic <- grepl("GS\\|", input$annotRow)
    am <- NULL
    if (any(ic)) {
      am <- cbind(am, vapply(input$annotRow[ic], function(i) {
        d <- gsdf()$featureId[gsdf()$gsId == sub("^GS\\|", 
                                                 "", i)]
        rownames(fd()) %fin% as.character(d)
      }, FUN.VALUE = logical(nrow(fd()))))
      colnames(am) <- input$annotRow[ic]
    }
    if (any(!ic)) {
      am2 <- fd()[, input$annotRow[!ic], drop = FALSE]
      if (is.null(am)) 
        am <- am2
      else am <- cbind(am, am2)
    }
    am <- as.data.frame(am[, input$annotRow, drop = FALSE])
    addHeatmapAnnotation(am[hm()$ord_r, , drop = FALSE], 
                         column = FALSE, var.name = input$annotRow)
  })
  output$rowSideCol <- renderPlot({
    par(mar = c(input$marginBottom, 0, 0, 0))
    addHeatmapAnnotation_plot(dat_rowSideCol(), ylim = ranges$y - 
                                0.5)
  })
  output$heatmap <- renderPlot({
    par(mar = c(input$marginBottom, 0, 0, input$marginRight))
    req(hm()$mat)
    image(hm()$mat, x = seq_len(nrow(hm()$mat)), y = seq_len(ncol(hm()$mat)), 
          xlim = ranges$x, ylim = ranges$y, col = heatColor(), 
          axes = FALSE, xlab = "", ylab = "", breaks = hm()$brk)
    if (ranges$x[2] - ranges$x[1] <= 60) {
      irn <- .rg(ranges$x, rownames(hm()$mat))
      mtext(side = 1, at = irn$at, text = irn$lab, las = 2, 
            line = 0.5)
      abline(v = seq(ranges$x[1], ranges$x[2], by = 1), 
             col = "white")
    }
    if (ranges$y[2] - ranges$y[1] <= 100) {
      rrn <- .rg(ranges$y, colnames(hm()$mat))
      abline(h = seq(ranges$y[1], ranges$y[2], by = 1), 
             col = "white")
      mtext(side = 4, at = rrn$at, text = rrn$lab, las = 2, 
            line = 0.5)
    }
  })
  output$dendCol <- renderPlot({
    req(hm()$dend_c)
    par(mar = c(0, 0, 1, input$marginRight))
    plot(hm()$dend_c, xaxs = "i", yaxs = "i", axes = FALSE, 
         xlim = ranges$x, center = TRUE)
    axis(side = 4)
  })
  output$dendRow <- renderPlot({
    req(hm()$dend_r)
    par(mar = c(input$marginBottom, 1, 0, 0))
    plot(hm()$dend_r, horiz = TRUE, yaxs = "i", xaxs = "i", 
         axes = FALSE, ylim = ranges$y)
    axis(side = 1)
  })
  ranges <- reactiveValues(x = NULL, y = NULL)
  observe({
    req(hm()$mat)
    if (!is.null(status()$ranges_x)) 
      ranges$x <- status()$ranges_x
    else ranges$x <- c(0, nrow(hm()$mat)) + 0.5
    if (!is.null(status()$ranges_y)) 
      ranges$y <- status()$ranges_y
    else ranges$y <- c(0, ncol(hm()$mat)) + 0.5
  })
  .rg <- function(x, tx) {
    v1 <- ceiling(x[1])
    v2 <- floor(x[2])
    list(at = v1:v2, lab = tx[v1:v2])
  }
  heatColor <- reactive({
    colorRampPalette(rev(brewer.pal(n = 7, name = input$heatmapColors)))(100)
  })
  observeEvent(input$heatmap_dblclick, {
    brush <- input$heatmap_brush
    if (!is.null(brush)) {
      ranges$x <- c(max(floor(brush$xmin) - 0.5, 0.5), 
                    min(round(brush$xmax) + 0.5, nrow(hm()$mat) + 
                          0.5))
      ranges$y <- c(max(floor(brush$ymin) - 0.5, 0.5), 
                    min(round(brush$ymax) + 0.5, ncol(hm()$mat) + 
                          0.5))
    }
    else {
      ranges$x <- c(0, nrow(hm()$mat)) + 0.5
      ranges$y <- c(0, ncol(hm()$mat)) + 0.5
    }
  })
  observeEvent(input$dendRow_dblclick, {
    brush <- input$dendRow_brush
    if (!is.null(brush)) {
      ranges$y <- c(max(floor(brush$ymin) - 0.5, 0.5), 
                    min(ceiling(brush$ymax) + 0.5, ncol(hm()$mat) + 
                          0.5))
    }
    else {
      ranges$y <- c(0, ncol(hm()$mat)) + 0.5
    }
  })
  observeEvent(input$dendCol_dblclick, {
    brush <- input$dendCol_brush
    if (!is.null(brush)) {
      ranges$x <- c(max(floor(brush$xmin) - 0.5, 0.5), 
                    min(ceiling(brush$xmax) + 0.5, nrow(hm()$mat) + 
                          0.5))
    }
    else {
      ranges$x <- c(0, nrow(hm()$mat)) + 0.5
    }
  })
  observeEvent(input$rowSideCol_dblclick, {
    brush <- input$rowSideCol_brush
    if (!is.null(brush)) {
      ranges$y <- c(max(floor(brush$ymin) - 0.5, 0.5), 
                    min(ceiling(brush$ymax) + 0.5, ncol(hm()$mat) + 
                          0.5))
    }
    else {
      ranges$y <- c(0, ncol(hm()$mat)) + 0.5
    }
  })
  observeEvent(input$colSideCol_dblclick, {
    brush <- input$colSideCol_brush
    if (!is.null(brush)) {
      ranges$x <- c(max(floor(brush$xmin) - 0.5, 0.5), 
                    min(ceiling(brush$xmax) + 0.5, nrow(hm()$mat) + 
                          0.5))
    }
    else {
      ranges$x <- c(0, nrow(hm()$mat)) + 0.5
    }
  })
  output$key_heatmap <- renderPlot(heatmapKey(range(hm()$mat, 
                                                    na.rm = TRUE), heatColor()))
  output$key_heatmap_ui <- renderUI({
    plotOutput(ns("key_heatmap"), height = "45px")
  })
  output$key_colSideCor <- renderPlot({
    lg <- dat_colSideCol()
    req(lg$key)
    if (lg$type == 2) {
      sideCorKey(x = lg$key, label = lg$var.name)
    }
    else {
      graphics::layout(matrix(seq_along(lg$key), 1))
      for (i in seq_along(lg$key)) sideCorKey(x = lg$key[[i]], 
                                              label = lg$var.name[i])
    }
  })
  output$key_colSideCor_ui <- renderUI({
    if (is.null(dat_colSideCol()$key)) 
      return()
    tagList(h4("Column bars"), plotOutput(ns("key_colSideCor"), 
                                          height = 266))
  })
  output$key_rowSideCor <- renderPlot({
    lg <- dat_rowSideCol()
    req(lg$key)
    if (lg$type == 2) {
      sideCorKey(x = lg$key, label = lg$var.name)
    }
    else {
      graphics::layout(matrix(seq_along(lg$key), 1))
      for (i in seq_along(lg$key)) sideCorKey(x = lg$key[[i]], 
                                              label = lg$var.name[i])
    }
  })
  output$key_rowSideCor_ui <- renderUI({
    if (is.null(dat_rowSideCol()$key)) 
      return()
    tagList(h4("Row bars"), plotOutput(ns("key_rowSideCor"), 
                                       height = 266))
  })
  output$empty1 <- renderPlot({
    par(mar = c(0, 0, 0, 0))
    plot(0, axes = FALSE, col = NA)
  })
  output$empty2 <- renderPlot({
    par(mar = c(0, 0, 0, 0))
    plot(0, axes = FALSE, col = NA)
  })
  output$empty3 <- renderPlot({
    par(mar = c(0, 0, 0, 0))
    plot(0, axes = FALSE, col = NA)
  })
  output$empty4 <- renderPlot({
    par(mar = c(0, 0, 0, 0))
    plot(0, axes = FALSE, col = NA)
  })
  show_col_dend <- reactiveVal(TRUE)
  observe(show_col_dend(input$colSortBy == "hierarchical cluster" || 
                          grepl("^HCL ", input$colSortBy)))
  show_col_sideColor <- reactiveVal(FALSE)
  observe(show_col_sideColor(length(input$annotCol) != 0))
  show_row_dend <- reactiveVal(TRUE)
  observe(show_row_dend(input$rowSortBy == "hierarchical cluster" || 
                          grepl("^HCL ", input$rowSortBy)))
  show_row_sideColor <- reactiveVal(FALSE)
  observe(show_row_sideColor(length(input$annotRow) != 0))
  width_left <- reactive(2)
  width_mid <- reactive({
    if (length(input$annotRow) <= 3) 
      r <- 1
    else if (length(input$annotRow) < 6) 
      r <- 2
    else r <- 3
  })
  width_right <- reactive({
    if (show_row_sideColor() & show_row_dend()) {
      wid <- 12 - width_mid() - width_left()
    }
    else if (!show_row_sideColor() & show_row_dend()) {
      wid <- 12 - width_left()
    }
    else if (show_row_sideColor() & !show_row_dend()) {
      wid <- 12 - width_mid()
    }
    else wid <- 12
    wid
  })
  height_colSideCol <- reactive({
    if (dat_colSideCol()$type == 1) 
      return("50px")
    paste0(20 * length(input$annotCol), "px")
  })
  output$topleft <- renderUI({
    if (!show_col_dend() || !show_row_dend()) 
      return(NULL)
    column(width_left(), plotOutput(ns("empty1"), height = "100px"), 
           offset = 0, style = "padding:0px;")
  })
  output$topmiddle <- renderUI({
    if (!show_col_dend() || !show_row_sideColor()) 
      return(NULL)
    column(width_mid(), plotOutput(ns("empty2"), height = "100px"), 
           offset = 0, style = "padding:0px;")
  })
  output$column_dend_ui <- renderUI({
    if (!show_col_dend()) 
      return(NULL)
    column(width_right(), plotOutput(ns("dendCol"), dblclick = ns("dendCol_dblclick"), 
                                     brush = brushOpts(id = ns("dendCol_brush"), resetOnNew = TRUE), 
                                     height = "100px"), offset = 0, style = "padding:0px;")
  })
  output$middleleft <- renderUI({
    if (!show_col_sideColor() || !show_row_dend()) 
      return(NULL)
    column(width_left(), plotOutput(ns("empty3"), height = height_colSideCol()), 
           offset = 0, style = "padding:0px;")
  })
  output$middlemiddle <- renderUI({
    if (!show_col_sideColor() || !show_row_sideColor()) 
      return(NULL)
    column(width_mid(), plotOutput(ns("empty4"), height = height_colSideCol()), 
           offset = 0, style = "padding:0px;")
  })
  output$column_sideColor_ui <- renderUI({
    if (!show_col_sideColor()) 
      return(NULL)
    column(width_right(), plotOutput(ns("colSideCol"), click = ns("colSideCol_click"), 
                                     dblclick = ns("colSideCol_dblclick"), hover = hoverOpts(id = ns("colSideCol_hover"), 
                                                                                             delay = 150), brush = brushOpts(id = ns("colSideCol_brush"), 
                                                                                                                             resetOnNew = TRUE), height = height_colSideCol()), 
           offset = 0, style = "padding:0px;")
  })
  output$row_dend_ui <- renderUI({
    if (!show_row_dend()) 
      return(NULL)
    column(width_left(), plotOutput(ns("dendRow"), dblclick = ns("dendRow_dblclick"), 
                                    brush = brushOpts(id = ns("dendRow_brush"), resetOnNew = TRUE), 
                                    height = "800px"), offset = 0, style = "padding:0px;")
  })
  output$row_sideColor_ui <- renderUI({
    if (!show_row_sideColor()) 
      return(NULL)
    column(width_mid(), plotOutput(ns("rowSideCol"), click = ns("rowSideCol_click"), 
                                   dblclick = ns("rowSideCol_dblclick"), hover = hoverOpts(id = ns("rowSideCol_hover"), 
                                                                                           delay = 150), brush = brushOpts(id = ns("rowSideCol_brush"), 
                                                                                                                           resetOnNew = TRUE), height = "800px"), offset = 0, 
           style = "padding:0px;")
  })
  output$heatmap_ui <- renderUI({
    column(width_right(), plotOutput(ns("heatmap"), click = ns("heatmap_click"), 
                                     dblclick = ns("heatmap_dblclick"), hover = hoverOpts(id = ns("heatmap_hover"), 
                                                                                          delay = 150), brush = brushOpts(id = ns("heatmap_brush"), 
                                                                                                                          resetOnNew = TRUE), height = "800px"), offset = 0, 
           style = "padding:0px;")
  })
  dat_tooltip <- reactive({
    req(pd())
    req(fd())
    list(pd = pd()[hm()$ord_c, colnames(pd()) %in% input$tooltipInfo, 
                   drop = FALSE], fd = fd()[hm()$ord_r, colnames(fd()) %in% 
                                              input$tooltipInfo, drop = FALSE])
  })
  tooltip_mouse <- reactiveVal(NULL)
  observeEvent(list(input$annotCol, input$annotRow), {
    tooltip_mouse(NULL)
  })
  observe({
    res <- NULL
    req(input$colSideCol_hover)
    req(input$colSideCol_click)
    hover_x <- ceiling(input$colSideCol_hover$x)
    hover_y <- ceiling(input$colSideCol_hover$y)
    x <- ceiling(input$colSideCol_click$x)
    y <- ceiling(input$colSideCol_click$y)
    if (hover_x == x && hover_y == y && dat_colSideCol()$type %in% 
        2:3) {
      if (dat_colSideCol()$type == 2) {
        varname <- dat_colSideCol()$var.name
        key <- dat_colSideCol()$key
        leg <- names(dat_colSideCol()$key$color)[x]
      }
      else {
        varname <- rownames(dat_colSideCol()$cmat)[y]
        key <- dat_colSideCol()$key[[match(varname, dat_colSideCol()$var.name)]]
        color <- match(dat_colSideCol()$cmat[y, x], dat_colSideCol()$cmat[y, 
        ])
        leg <- names(key$color)[color]
      }
      res <- sprintf("<b>%s:</b> %s", varname, leg)
    }
    tooltip_mouse(res)
  })
  observe({
    res <- NULL
    req(input$rowSideCol_hover)
    req(input$rowSideCol_click)
    hover_x <- ceiling(input$rowSideCol_hover$x)
    hover_y <- ceiling(input$rowSideCol_hover$y)
    x <- ceiling(input$rowSideCol_click$x)
    y <- ceiling(input$rowSideCol_click$y)
    if (hover_x == x && hover_y == y && dat_rowSideCol()$type %in% 
        2:3) {
      if (dat_rowSideCol()$type == 2) {
        varname <- dat_rowSideCol()$var.name
        key <- dat_rowSideCol()$key
        leg <- names(dat_rowSideCol()$key$color)[y]
      }
      else {
        varname <- rownames(dat_rowSideCol()$cmat)[x]
        key <- dat_rowSideCol()$key[[match(varname, dat_rowSideCol()$var.name)]]
        color <- match(dat_rowSideCol()$cmat[x, y], dat_rowSideCol()$cmat[x, 
        ])
        leg <- names(key$color)[color]
      }
      res <- sprintf("<b>%s:</b> %s", varname, leg)
    }
    tooltip_mouse(res)
  })
  clickedName <- reactiveVal(NULL)
  observe({
    res <- NULL
    req(input$heatmap_hover)
    req(input$heatmap_click)
    hover_x <- round(input$heatmap_hover$x)
    hover_y <- round(input$heatmap_hover$y)
    x <- round(input$heatmap_click$x)
    y <- round(input$heatmap_click$y)
    if (hover_x == x && hover_y == y) {
      l <- c(dat_tooltip()$pd[x, , drop = FALSE], dat_tooltip()$fd[y, 
                                                                   , drop = FALSE])
      tl <- vapply(names(l), function(x) {
        if (is.numeric(cv <- l[[x]])) 
          cv <- round(cv, digits = 2)
        sprintf("<b>%s:</b> %s", x, cv)
      }, character(1))
      tl <- paste(tl, collapse = "<br>")
      res <- sprintf("<b>Sample:</b> %s <br> <b>Feature:</b> %s", 
                     rownames(hm()$mat)[x], colnames(hm()$mat)[y])
      res <- paste(res, tl, sep = "<br>")
    }
    clickedName(c(col = rownames(hm()$mat)[x], row = colnames(hm()$mat)[y]))
    tooltip_mouse(res)
  })
  callModule(shinyPlotTooltips, id = "tlp_heatmap", points = tooltip_mouse)
  brushedValues <- reactive({
    brush <- input$heatmap_brush
    l <- list(col = NULL, row = NULL)
    if (!is.null(brush)) {
      x1 <- ceiling(max(round(brush$xmin), 0.5))
      x2 <- floor(min(round(brush$xmax), nrow(hm()$mat)))
      y1 <- ceiling(max(round(brush$ymin), 0.5))
      y2 <- floor(min(round(brush$ymax), ncol(hm()$mat)))
      l <- list(row = colnames(hm()$mat)[y1:y2], col = rownames(hm()$mat)[x1:x2])
    }
    l
  })
  selVal <- reactiveValues(clicked = NULL, selected = list(col = NULL, 
                                                           row = NULL))
  observeEvent(input$clear, {
    selVal$clicked <- NULL
    selVal$selected <- list(col = NULL, row = NULL)
  })
  observeEvent(list(clickedName(), brushedValues()), {
    selVal$clicked <- clickedName()
    selVal$selected <- brushedValues()
  })
  reactive({
    r <- list(clicked = selVal$clicked, brushed = selVal$selected)
    attr(r, "status") <- list(annotCol = input$annotCol, 
                              annotRow = input$annotRow, colSortBy = input$colSortBy, 
                              rowSortBy = input$rowSortBy, tooltipInfo = input$tooltipInfo, 
                              marginRight = input$marginRight, marginBottom = input$marginBottom, 
                              heatmapColors = input$heatmapColors, scale = input$scale, 
                              clusterColDist = input$clusterColDist, clusterColLink = input$clusterColLink, 
                              clusterRowDist = input$clusterRowDist, clusterRowLink = input$clusterRowLink, 
                              rowDendrogram = rowSB()$hcl, rowOrder = rowSB()$ord, 
                              colDendrogram = colSB()$hcl, colOrder = colSB()$ord, 
                              ranges_x = ranges$x, ranges_y = ranges$y)
    r
  })
}


iheatmapOutput <- function (id) 
{
  ns <- NS(id)
  tagList(fluidRow(uiOutput(ns("topleft")), uiOutput(ns("topmiddle")), 
                   uiOutput(ns("column_dend_ui")), uiOutput(ns("middleleft")), 
                   uiOutput(ns("middlemiddle")), uiOutput(ns("column_sideColor_ui")), 
                   uiOutput(ns("row_dend_ui")), uiOutput(ns("row_sideColor_ui")), 
                   shinycssloaders::withSpinner(uiOutput(ns("heatmap_ui")), 
                                                type = 8, color = "green")), shinyPlotTooltipsUI(ns("tlp_heatmap")))
}


L1_data_space_module <- function (input, output, session, expr, pdata, fdata, reactive_x_s = reactive(NULL), 
                                  reactive_y_s = reactive(NULL), reactive_x_f = reactive(NULL), 
                                  reactive_y_f = reactive(NULL), cormat = reactive(NULL), status = reactive(NULL), imputationMethod = 'custom') 
{
  ns <- session$ns
  cmat <- reactive({
    req(expr())
    if (ncol(expr()) <= 3) 
      return(NULL)
    if (is.null(cormat())) {
      cc <- cor(expr(), use = "pairwise.complete.obs")
      diag(cc) <- NA
      hcl <- hclust(as.dist(1 - cor(t(cc), use = "pair")), 
                    method = "ward.D")
      dend <- as.dendrogram(hcl)
      dl <- list(pearson_ward.D = list(ord = hcl$order, 
                                       hcl = dend))
      attr(cc, "rowDendrogram") <- dl
      attr(cc, "colDendrogram") <- dl
      return(cc)
    }
    else if (nrow(cormat()) == ncol(cormat()) && nrow(cormat()) == 
             ncol(expr())) {
      return(cormat())
    }
    else stop("Incorrect dimension of cormat!")
  })
  
  s_cor_heatmap <- callModule(iheatmapModule, "corheatmapViewer", 
                              mat = cmat,
                              pd = pdata,
                              fd = pdata,
                              status = reactive(status()$eset_cor_heatmap), 
                              fill.NA = FALSE
                              )

  s_heatmap <- callModule(iheatmapModule,
                           "heatmapViewer", 
                            mat = expr,
                            pd = pdata,
                            fd = fdata,
                            status = reactive(status()$eset_heatmap),
                            fill.NA = TRUE,
                            method = imputationMethod
                            )

  s_feature_fig <- callModule(meta_scatter_module, id = "feature_space", 
                              reactive_meta = fdata, reactive_expr = expr, combine = "feature", 
                              source = "scatter_meta_feature", reactive_x = reactive_x_f, 
                              reactive_y = reactive_y_f, reactive_status = reactive(status()$eset_fdata_fig))
  s_sample_fig <- callModule(meta_scatter_module, id = "sample_space", 
                             reactive_meta = pdata, reactive_expr = expr, combine = "pheno", 
                             source = "scatter_meta_sample", reactive_x = reactive_x_s, 
                             reactive_y = reactive_y_s, reactive_status = reactive(status()$eset_pdata_fig))
  tab_rows_fdata <- reactiveVal(TRUE)
  tab_rows_pdata <- reactiveVal(TRUE)
  notNullAndPosLength <- function(x) !is.null(x) && length(x) > 
    0
  observeEvent(s_cor_heatmap(), {
    if (notNullAndPosLength(s_cor_heatmap()$brushed$col)) {
      tab_rows_pdata(s_cor_heatmap()$brushed$col)
    }
    else if (notNullAndPosLength(s_cor_heatmap()$clicked)) {
      tab_rows_pdata(s_cor_heatmap()$clicked[["col"]])
    }
    else tab_rows_pdata(TRUE)
  })
  observeEvent(s_heatmap(), {
    if (notNullAndPosLength(s_heatmap()$brushed$row)) {
      tab_rows_fdata(s_heatmap()$brushed$row)
    }
    else if (notNullAndPosLength(s_heatmap()$clicked)) {
      tab_rows_fdata(s_heatmap()$clicked[["row"]])
    }
    else {
      tab_rows_fdata(TRUE)
    }
    if (notNullAndPosLength(s_heatmap()$brushed$col)) {
      tab_rows_pdata(s_heatmap()$brushed$col)
    }
    else if (notNullAndPosLength(s_heatmap()$clicked)) {
      tab_rows_pdata(s_heatmap()$clicked[["col"]])
    }
    else {
      tab_rows_pdata(TRUE)
    }
  })
  observeEvent(c(s_feature_fig()), {
    if (notNullAndPosLength(s_feature_fig()$selected)) {
      tab_rows_fdata(s_feature_fig()$selected)
    }
    else if (notNullAndPosLength(s_feature_fig()$clicked)) {
      tab_rows_fdata(s_feature_fig()$clicked)
    }
    else tab_rows_fdata(TRUE)
  })
  observeEvent(s_sample_fig(), {
    if (notNullAndPosLength(s_sample_fig()$selected)) {
      tab_rows_pdata(s_sample_fig()$selected)
    }
    else if (notNullAndPosLength(s_sample_fig()$clicked)) {
      tab_rows_pdata(s_sample_fig()$clicked)
    }
    else tab_rows_pdata(TRUE)
  })
  tab_pd <- callModule(dataTable_module, id = "tab_pheno", 
                       reactive_data = pdata, tab_status = reactive(status()$eset_pdata_tab), 
                       tab_rows = tab_rows_pdata)
  tab_fd <- callModule(dataTable_module, id = "tab_feature", 
                       reactive_data = fdata, tab_status = reactive(status()$eset_fdata_tab), 
                       tab_rows = tab_rows_fdata)
  tab_expr <- callModule(dataTable_module, id = "tab_expr", 
                         reactive_data = reactive({
                           req(expr())
                           cbind(data.frame(feature = rownames(expr()), expr()))
                         }), tab_status = reactive(status()$eset_exprs_tab), tab_rows = tab_rows_fdata, 
                         selector = FALSE)
  selectedFeatures <- reactiveVal()
  selectedSamples <- reactiveVal()
  observeEvent(s_cor_heatmap(), {
    if (!is.null(s_cor_heatmap()$brushed$col)) {
      selectedSamples(s_cor_heatmap()$brushed$col)
    }
    else if (!is.null(s_cor_heatmap()$clicked)) {
      selectedSamples(s_cor_heatmap()$clicked["col"])
    }
    else selectedSamples(character(0))
  })
  observeEvent(s_heatmap(), {
    if (!is.null(s_heatmap()$brushed$row)) {
      selectedFeatures(s_heatmap()$brushed$row)
    }
    else if (!is.null(s_heatmap()$clicked)) {
      selectedFeatures(s_heatmap()$clicked["row"])
    }
    else selectedFeatures(character(0))
    if (!is.null(s_heatmap()$brushed$col)) {
      selectedSamples(s_heatmap()$brushed$col)
    }
    else if (!is.null(s_heatmap()$clicked)) {
      selectedSamples(s_heatmap()$clicked["col"])
    }
    else selectedSamples(character(0))
  })
  observeEvent(s_feature_fig(), {
    if (!is.null(s_feature_fig()$selected) && length(s_feature_fig()$selected) > 
        0) {
      selectedFeatures(s_feature_fig()$selected)
    }
    else if (!is.null(s_feature_fig()$clicked)) 
      selectedFeatures(s_feature_fig()$clicked)
  })
  observeEvent(s_sample_fig(), {
    if (!is.null(s_sample_fig()$selected) && length(s_sample_fig()$selected) > 
        0) {
      selectedSamples(s_sample_fig()$selected)
    }
    else if (!is.null(s_sample_fig()$clicked)) 
      selectedSamples(s_sample_fig()$clicked)
  })
  observeEvent(tab_pd(), {
    selectedSamples(tab_pd())
  })
  singleTrue <- function(x) !is.null(x) && is.logical(x) && 
    length(x) == 1 && x
  observeEvent(tab_fd(), {
    i1 <- length(tab_fd()) < length(tab_rows_fdata())
    i2 <- singleTrue(tab_rows_fdata())
    i3 <- !singleTrue(tab_fd())
    if ((i1 || i2) && i3) 
      selectedFeatures(tab_fd())
  })
  observeEvent(tab_expr(), {
    i1 <- length(tab_expr()) < length(tab_rows_fdata())
    i2 <- singleTrue(tab_rows_fdata())
    i3 <- !singleTrue(tab_expr())
    if ((i1 || i2) && i3) 
      selectedFeatures(tab_expr())
  })
  tab_gslist <- callModule(gslist_module, id = "gsList", reactive_i = tab_rows_fdata, 
                           reactive_featureData = fdata)
  observeEvent(tab_gslist(), {
    req(tab_gslist())
    selectedFeatures(tab_gslist())
  })
  observe({
    if (!is.null(tb <- status()$eset_active_tab)) 
      updateNavbarPage(session = session, inputId = "eset", 
                       selected = tb)
  })
  na2null <- function(x) {

    if (is.null(x) || is.na(x) || length(x) == 0) 
      return(NULL)
    x
  }
  observe({
    tab_rows_fdata(status()$eset_fdata_tabrows)
    tab_rows_pdata(status()$eset_pdata_tabrows)
    selectedSamples(na2null(status()$eset_selected_samples))
    selectedFeatures(status()$eset_selected_features)
  })
  reactive({
    l <- list(feature = selectedFeatures(), sample = selectedSamples())
    sta <- list(eset_active_tab = input$eset, eset_pdata_tab = attr(tab_pd(), 
                                                                    "status"), eset_fdata_tab = attr(tab_fd(), "status"), 
                eset_exprs_tab = attr(tab_expr(), "status"), eset_fdata_fig = attr(s_feature_fig(), 
                                                                                   "status"), eset_pdata_fig = attr(s_sample_fig(), 
                                                                                                                    "status"), eset_heatmap = attr(s_heatmap(), "status"), 
                eset_fdata_tabrows = tab_rows_fdata(), eset_pdata_tabrows = tab_rows_pdata(), 
                eset_selected_samples = c(selectedSamples()), eset_selected_features = c(selectedFeatures()))
    if (sta$eset_active_tab != "Feature table") 
      sta$eset_fdata_tab$rows_selected <- NULL
    if (sta$eset_active_tab != "Sample table") 
      sta$eset_pdata_tab$rows_selected <- NULL
    attr(l, "status") <- sta
    l
  })
}


L1_data_space_ui <- function (id, activeTab = "Feature") 
{
  ns <- NS(id)
  navbarPage("Data", id = ns("eset"), selected = activeTab, 
             theme = shinytheme("spacelab"), tabPanel("Feature", meta_scatter_ui(ns("feature_space"))), 
             tabPanel("Feature table", dataTable_ui(ns("tab_feature"))), 
             tabPanel("Sample", meta_scatter_ui(ns("sample_space"))), 
             tabPanel("Sample table", dataTable_ui(ns("tab_pheno"))), 
             tabPanel("Cor", fluidRow(column(6, dropdown(inputId = "mydropdown2", 
                                                         label = "Controls", circle = FALSE, status = "default", 
                                                         icon = icon("cog"), width = 700, tooltip = tooltipOptions(title = "Click to update heatmap and check legend!"), 
                                                         margin = "10px", tabsetPanel(tabPanel("Parameters", 
                                                                                               iheatmapInput(id = ns("corheatmapViewer"), scaleOn = "none")), 
                                                                                      tabPanel("Legend", iheatmapLegend(id = ns("corheatmapViewer")))))), 
                                      column(6, align = "right", iheatmapClear(id = ns("corheatmapViewer"))), 
                                      column(12, iheatmapOutput(id = ns("corheatmapViewer"))))), 
             tabPanel("Heatmap", fluidRow(column(6, dropdown(inputId = "mydropdown", 
                                                             label = "Controls", circle = FALSE, status = "default", 
                                                             icon = icon("cog"), width = 700, tooltip = tooltipOptions(title = "Click to update heatmap and check legend!"), 
                                                             margin = "10px", tabsetPanel(tabPanel("Parameters", 
                                                                                                   iheatmapInput(id = ns("heatmapViewer"))), tabPanel("Legend", 
                                                                                                                                                      iheatmapLegend(id = ns("heatmapViewer")))))), 
                                          column(6, align = "right", iheatmapClear(id = ns("heatmapViewer"))), 
                                          column(12, iheatmapOutput(id = ns("heatmapViewer"))))), 
             tabPanel("Expression", dataTable_ui(ns("tab_expr"))), 
             tabPanel("GSList", gslist_ui(ns("gsList"))))
}


L1_result_space_module <- function (input, output, session, reactive_expr, reactive_phenoData, 
                                    reactive_featureData, reactive_i = reactive(NULL), reactive_highlight = reactive(NULL), 
                                    additionalTabs = NULL, object = NULL, status = reactive(NULL)) 
{
  ns <- session$ns
  v <- callModule(feature_general_module, id = "feature_general", 
                  reactive_expr = reactive_expr, reactive_i = reactive_i, 
                  reactive_highlight = reactive_highlight, reactive_phenoData = reactive_phenoData, 
                  reactive_featureData = reactive_featureData, reactive_status = reactive(status()$analyst_feature_general))
  v2 <- callModule(enrichment_fgsea_module, id = "fgsea", reactive_featureData = reactive_featureData, 
                   reactive_status = reactive(status()$analyst_fgsea))
  v3 <- callModule(enrichment_analysis_module, id = "ora", 
                   reactive_i = reactive_i, reactive_featureData = reactive_featureData)
  v4 <- callModule(string_module, id = "stringdb", reactive_ids = reactive({
    i <- grep("^StringDB\\|", colnames(reactive_featureData()))
    reactive_featureData()[reactive_i(), i[1]]
  }), reactive_status = reactive(status()$analyst_stringdb), 
  active = reactive(status()$analyst_active_tab == "StringDB"))
  v5 <- callModule(sample_general_module, id = "sample_general", 
                   reactive_phenoData = reactive_phenoData, reactive_j = reactive_highlight, 
                   reactive_status = reactive(status()$analyst_sample_general))
  v6 <- callModule(geneshot_module, id = "geneshotTab", fdata = reactive_featureData, 
                   feature_selected = reactive_i, reactive_status = reactive(status()$analyst_gene_shot))
  v7 <- callModule(ptmotif_module, id = "ptm", fdata = reactive_featureData, 
                   feature_selected = reactive_i)
  if (length(additionalTabs) > 0) {
    for (lo in additionalTabs) {
      callModule(lo$moduleServer, id = lo$moduleName, pdata = reactive_phenoData, 
                 fdata = reactive_featureData, expr = reactive_expr, 
                 feature_selected = reactive_i, sample_selected = reactive_highlight, 
                 object = object)
    }
  }
  observe({
    if (!is.null(tb <- status()$analyst_active_tab)) 
      updateNavbarPage(session = session, inputId = "analyst", 
                       selected = tb)
  })
  output$optTabs <- renderUI({
    titleTabs <- list(title = "Analysis", id = ns("analyst"), 
                      theme = shinytheme("spacelab"), tabPanel("Feature", 
                                                               feature_general_ui(ns("feature_general"))))
    sampleAnalyst <- list(tabPanel("Sample", sample_general_ui(ns("sample_general"))))
    geneshot <- list(tabPanel("Geneshot", geneshot_ui(ns("geneshotTab"))))
    optionalTabs <- list()
    if (!is.null(attr(reactive_featureData(), "GS"))) {
      optionalTabs <- c(optionalTabs, list(tabPanel("ORA", 
                                                    enrichment_analysis_ui(ns("ora")))))
      optionalTabs <- c(optionalTabs, list(tabPanel("fGSEA", 
                                                    enrichment_fgsea_ui(ns("fgsea")))))
    }
    if (any(grepl("^StringDB\\|", colnames(reactive_featureData())))) 
      optionalTabs <- c(optionalTabs, list(tabPanel("StringDB", 
                                                    string_ui(ns("stringdb")))))
    if (any(grepl("^SeqLogo\\|", colnames(reactive_featureData())))) 
      optionalTabs <- c(optionalTabs, list(tabPanel("SeqLogo", 
                                                    ptmotif_ui(ns("ptm")))))
    if (length(additionalTabs) > 0) {
      for (lo in additionalTabs) {
        optionalTabs <- c(optionalTabs, list(tabPanel(lo$tabName, 
                                                      lo$moduleUi(ns(lo$moduleName)))))
      }
    }
    do.call(navbarPage, c(titleTabs, optionalTabs, geneshot, 
                          sampleAnalyst))
  })
  reactive({
    list(analyst_active_tab = input$analyst, analyst_feature_general = v(), 
         analyst_sample_general = v5(), analyst_gene_shot = v6(), 
         analyst_fgsea = v2(), analyst_stringdb = v4())
  })
}


L1_result_space_ui <- function (id) 
{
  ns <- NS(id)
  tagList(uiOutput(ns("optTabs")))
}


line_rect <- function (l, coord) 
{
  x <- l$x
  y <- l$y
  if (is.null(x) && is.null(y)) 
    return(NULL)
  if (!is.numeric(x) && !is.numeric(y)) 
    return(NULL)
  if (!is.numeric(x) && is.numeric(y) && !l$corner %in% c("top", 
                                                          "bottom")) 
    return(NULL)
  if (is.numeric(x) && !is.numeric(y) && !l$corner %in% c("left", 
                                                          "right")) 
    return(NULL)
  if (!is.numeric(coord$x) && is.numeric(coord$y) && !l$corner %in% 
      c("top", "bottom")) 
    return(NULL)
  if (is.numeric(coord$x) && !is.numeric(coord$y) && !l$corner %in% 
      c("left", "right")) 
    return(NULL)
  if (l$corner == "None") 
    return(NULL)
  minx <- min(coord$x, na.rm = TRUE)
  maxx <- max(coord$x, na.rm = TRUE)
  miny <- min(coord$y, na.rm = TRUE)
  maxy <- max(max(coord$y, na.rm = TRUE), y)
  x.exp <- (maxx - minx) * 0.05
  y.exp <- (maxy - miny) * 0.05
  if (x.exp == 0 || y.exp == 0) 
    return(NULL)
  rect <- NULL
  if (l$corner == "volcano") {
    ax <- abs(l$x)
    y0 <- max(l$y, miny - y.exp)
    y1 <- maxy + y.exp
    l_x0 <- minx - x.exp
    l_x1 <- min(-ax, maxx + x.exp)
    r_x0 <- max(ax, minx - x.exp)
    r_x1 <- maxx + x.exp
    rect <- list(c(x0 = l_x0, x1 = l_x1, y0 = y0, y1 = y1), 
                 c(x0 = r_x0, x1 = r_x1, y0 = y0, y1 = y1))
  }
  else {
    if (grepl("top", l$corner)) {
      y0 <- max(l$y, miny - y.exp)
      y1 <- maxy + y.exp
    }
    else if (grepl("bottom", l$corner)) {
      y0 <- miny - y.exp
      y1 <- min(l$y, maxy + y.exp)
    }
    else {
      y0 <- miny - y.exp
      y1 <- maxy + y.exp
    }
    if (grepl("left", l$corner)) {
      x0 <- minx - x.exp
      x1 <- min(l$x, maxx + x.exp)
    }
    else if (grepl("right", l$corner)) {
      x0 <- max(l$x, minx - x.exp)
      x1 <- maxx + x.exp
    }
    else {
      x0 <- minx - x.exp
      x1 <- maxx + x.exp
    }
    if ((x0 > x1) || (y0 > y1)) 
      x <- y <- x0 <- x1 <- y0 <- y1 <- NULL
    rect <- list(c(x0 = x0, x1 = x1, y0 = y0, y1 = y1))
  }
  list(x = x, y = y, rect = rect)
}


list2csc <- function (l, dimnames) 
{
  if (!is.list(dimnames)) 
    stop("list2csc: dimnames should be a list!")
  if (is.null(dimnames[[1]])) 
    stop("list2csc: dimnames should contain as least one element for the rownames of the matrix")
  if (!all(colnames(l) %in% c("featureId", "gsId", "weight"))) 
    stop("list2csc: colnames of l should be \"featureId\", \"gsId\", \"weight\"")
  names <- dimnames[[1]]
  if (length(dimnames) > 1) 
    cn <- dimnames[[2]]
  else cn <- unique(l$gsId)
  args <- list(i = c(fmatch(l$featureId, names)), j = c(fmatch(l$gsId, 
                                                               cn)), dims = c(length(names), length(cn)), dimnames = list(names, 
                                                                                                                          cn))
  if ("weight" %in% colnames(l)) 
    args$x <- l$weight
  do.call(sparseMatrix, args)
}


meta_scatter_module <- function (input, output, session, reactive_meta = reactive(NULL), 
                                 reactive_expr = reactive(NULL), combine = c("pheno", "feature"), 
                                 source = "plotlyscattersource", reactive_x = reactive(NULL), 
                                 reactive_y = reactive(NULL), reactive_status = reactive(NULL)) 
{
  ns <- session$ns
  triset <- reactive({
    ts <- trisetter(expr = reactive_expr(), meta = reactive_meta(), 
                    combine = combine[1])
    ts[ts[, 1] != "Surv", ]
  })
  xax <- reactiveVal()
  observeEvent(reactive_x(), {
    r <- list()
    if (!is.null(reactive_x())) {
      l <- strsplit(reactive_x(), "\\|")[[1]]
      r <- list(v1 = l[1], v2 = l[2], v3 = l[3])
    }
    xax(r)
  })
  yax <- reactiveVal()
  observe({
    r <- list()
    if (!is.null(reactive_y())) {
      l <- strsplit(reactive_y(), "\\|")[[1]]
      r <- list(v1 = l[1], v2 = l[2], v3 = l[3])
    }
    yax(r)
  })
  v1 <- callModule(triselector_module, id = "tris_main_scatter1", 
                   reactive_x = triset, label = "X-axis", reactive_selector1 = reactive(xax()$v1), 
                   reactive_selector2 = reactive(xax()$v2), reactive_selector3 = reactive(xax()$v3))
  v2 <- callModule(triselector_module, id = "tris_main_scatter2", 
                   reactive_x = triset, label = "Y-axis", reactive_selector1 = reactive(yax()$v1), 
                   reactive_selector2 = reactive(yax()$v2), reactive_selector3 = reactive(yax()$v3))
  pre_vol <- reactiveVal(FALSE)
  observe({
    vv <- c("v1", "v2", "v3")
    if (all(vv %in% names(xax())) && all(vv %in% names(yax()))) {
      if (xax()$v1 == "ttest" && yax()$v1 == "ttest" && 
          xax()$v3 == "mean.diff" && yax()$v3 %in% c("log.fdr", 
                                                     "log.pvalue")) 
        pre_vol(TRUE)
    }
  })
  attr4select_status <- reactiveVal()
  attr4select <- callModule(attr4selector_module, id = "a4selector", 
                            reactive_meta = reactive_meta, reactive_expr = reactive_expr, 
                            reactive_triset = triset, pre_volcano = pre_vol, reactive_status = attr4select_status)
  xycoord <- reactive({
    req(v1()$variable)
    req(v2()$variable)
    req(!v1()$variable %in% c("Select a variable!", ""))
    req(!v2()$variable %in% c("Select a variable!", ""))
    x <- varSelector(v1(), reactive_expr(), reactive_meta())
    y <- varSelector(v2(), reactive_expr(), reactive_meta())
    req(x)
    req(y)
    req(is.numeric(x) || is.numeric(y))
    req(length(x) == length(y))
    list(x = x, y = y)
  })
  rectval <- reactiveVal(NULL)
  observe({
    req(xycoord())
    rectval(line_rect(l = attr4select$cutoff, xycoord())$rect)
  })
  observeEvent(input$clear, {
    pre_vol(FALSE)
    rectval(NULL)
  })
  observeEvent(list(v1(), v2()), {
    if (is.null(attr4select$cutoff) || attr4select$cutoff$corner == 
        "None") 
      rectval(NULL)
  })
  scatter_vars <- reactive({
    req(l <- xycoord())
    l$source <- source
    l$xlab <- attr(l$x, "label")
    l$ylab <- attr(l$y, "label")
    l$color <- attr4select$color
    l$shape <- attr4select$shape
    l$size <- attr4select$size
    l$tooltips <- attr4select$tooltips
    l$highlight <- attr4select$highlight
    l$highlightName <- attr4select$highlightName
    l$rect <- rectval()
    l$inSelection <- NA
    l
  })
  showRegLine <- reactiveVal(FALSE)
  htestV1 <- reactiveVal()
  htestV2 <- reactiveVal()
  v_scatter <- callModule(plotly_scatter_module, id = "main_scatterOutput", 
                          reactive_param_plotly_scatter = scatter_vars, reactive_regLine = showRegLine, 
                          htest_var1 = htestV1, htest_var2 = htestV2)
  observe({
    showRegLine(v_scatter()$regline)
  })
  selVal <- reactiveVal(list(clicked = character(0), selected = character(0)))
  sbc <- reactiveVal(FALSE)
  observeEvent(list(input$clear, reactive_expr()), {
    selVal(list(clicked = character(0), selected = character(0)))
    sbc(FALSE)
  })
  clientSideSelection <- reactiveVal(character(0))
  observeEvent(v_scatter(), {
    if (combine == "pheno") 
      l <- colnames(reactive_expr())
    else l <- rownames(reactive_expr())
    u_c <- l[v_scatter()$clicked]
    u_s <- l[v_scatter()$selected]
    req(!identical(tmp <- c(u_c, u_s), clientSideSelection()))
    clientSideSelection(tmp)
    selVal(list(clicked = u_c, selected = u_s))
    sbc(FALSE)
  })
  emptyValues <- function(...) {
    l <- list(...)
    j <- vapply(l, function(x) is.null(x) || length(x) == 
                  0, logical(1))
    all(j)
  }
  returnCornerSelection <- reactiveVal(TRUE)
  observeEvent(rectval(), {
    if (!returnCornerSelection()) 
      return(NULL)
    rec <- rectval()
    if (is.null(rec)) {
      selVal(list(clicked = character(0), selected = character(0)))
      return(NULL)
    }
    req(cc <- xycoord())
    if (combine == "pheno") 
      l <- colnames(reactive_expr())
    else l <- rownames(reactive_expr())
    i <- lapply(rec, function(r1) {
      which(cc$x > r1["x0"] & cc$x < r1["x1"] & cc$y > 
              r1["y0"] & cc$y < r1["y1"])
    })
    i <- sort(unique(unlist(i)))
    selVal(list(clicked = character(0), selected = l[i]))
    sbc(TRUE)
  })
  observe({
    sv <- selVal()
    attr(sv, "status") <- list(xax = v1(), yax = v2(), showRegLine = showRegLine(), 
                               attr4 = attr4select$status, htestV1 = v_scatter()$htest_V1, 
                               htestV2 = v_scatter()$htest_V2, selection_clicked = sv$clicked, 
                               selection_selected = sv$selected, selectByCorner = sbc())
    selVal(sv)
  })
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status())) 
      return()
    xax(NULL)
    xax(list(v1 = s$xax[[1]], v2 = s$xax[[2]], v3 = s$xax[[3]]))
  })
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status())) 
      return()
    yax(NULL)
    yax(list(v1 = s$yax[[1]], v2 = s$yax[[2]], v3 = s$yax[[3]]))
  })
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status())) 
      return()
    attr4select_status(NULL)
    attr4select_status(s$attr4)
  })
  observeEvent(reactive_status(), {
    if (!is.null(s <- reactive_status())) 
      showRegLine(s$showRegLine)
  })
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status())) 
      return()
    htestV1(s$htestV1)
    htestV2(s$htestV2)
  })
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status())) 
      return()
    returnCornerSelection(s$selectByCorner)
    selVal(list(clicked = s$selection_clicked, selected = s$selection_selected))
  })
  selVal
}


meta_scatter_ui <- function (id) 
{
  ns <- NS(id)
  tagList(fluidRow(column(1, attr4selector_ui(ns("a4selector")), 
                          actionBttn(ns("clear"), "Clear figure selection", style = "minimal", 
                                     color = "primary", size = "xs")), column(11, triselector_ui(ns("tris_main_scatter1")), 
                                                                              triselector_ui(ns("tris_main_scatter2")))), plotly_scatter_ui(ns("main_scatterOutput"), 
                                                                                                                                            height = "666px"))
}


motifRF <- function (fg.seqs, bg.seqs, fg.pfm = NULL, bg.pfm = NULL) 
{
  aa <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
          "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  if (!is.null(fg.pfm)) 
    fg <- fg.pfm
  else fg <- aaFreq(fg.seqs)
  if (!is.null(bg.pfm)) 
    bk <- bg.pfm
  else bk <- aaFreq(bg.seqs)
  motif <- fg/bk
  motif[is.na(motif)] <- 0
  motif <- sweep(motif, 2, colSums(motif), "/")
  rownames(motif) <- aa
  motif
}


multi.t.test <- function (x, pheno, compare = NULL, fillNA = FALSE, method = 'perseus', ...) 
{
  print(paste0(' #### using the imputation method from multi.t.test #### ',method))
  x0 <- x

  if (is.vector(compare) || length(compare) == 3) 
    compare <- matrix(compare, nrow = 1)
  if (fillNA) 
    x <- fillNA(x, method = method)
  if (is.null(compare)) 
    return(NULL)
  tl <- lapply(unique(compare[, 1]), function(x) {
    x <- compare[x == compare[, 1], -1, drop = FALSE]
    unique(unlist(split(x, row(x))))
  })
  names(tl) <- unique(compare[, 1])
  df <- data.frame(row.names = rownames(x))
  for (i in names(tl)) {
    for (j in tl[[i]]) {
      m <- x[, which(pheno[[i]] == j), drop = FALSE]
      m0 <- x0[, which(pheno[[i]] == j), drop = FALSE]
      rv <- rowSums(!is.na(m0))
      rm <- rowMeans(m, na.rm = TRUE)
      rq <- rank(rm, na.last = TRUE)/sum(!is.na(rm))
      rq[is.na(rm)] <- NA
      df[[paste("mean", i, j, sep = "|")]] <- rm
      df[[paste("n value", i, j, sep = "|")]] <- rv
      df[[paste("quantile", i, j, sep = "|")]] <- rq
    }
  }
  for (i in seq_len(nrow(compare))) {
    v <- compare[i, ]
    i1 <- which(pheno[[v[1]]] == v[2])
    i2 <- which(pheno[[v[1]]] == v[3])
    if (length(i1) == 0 || length(i2) == 0) 
      stop("Didn't find var: ", paste(v, collapse = "-"))
    tv <- apply(x, 1, function(xx) {
      t <- try(t.test(xx[i1], xx[i2], var.equal = TRUE, 
                      ...), silent = TRUE)
      if (!is(t, "htest")) 
        return(c(pvalue = NA, mean.diff = NA))
      if (length(t$estimate) == 1) 
        md <- t$estimate[[1]]
      else md <- t$estimate[1] - t$estimate[2]
      c(pvalue = t$p.value, mean.diff = md)
    })
    pv <- tv[1, ]
    fdr <- p.adjust(pv, method = "fdr")
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "pvalue", 
              sep = "|")]] <- pv
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "log.pvalue", 
              sep = "|")]] <- -log10(pv)
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "fdr", 
              sep = "|")]] <- fdr
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "log.fdr", 
              sep = "|")]] <- -log10(fdr)
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "mean.diff", 
              sep = "|")]] <- tv[2, ]
  }
  df
}


nColors <- function (k, stop = FALSE) 
{
  k <- as.integer(k)
  if (!is.integer(k) || k > 60 || k < 1 || is.na(k)) {
    if (stop) 
      stop("k should be an integer between 1 and 60!")
    else return(NULL)
  }
  l <- list(c("#B5D3A2"), c("#C07CCA", "#B5D3A2"), c("#B0C8C8", 
                                                     "#C470C8", "#B6D778"), c("#9FD2D2", "#D8949E", "#B86BD9", 
                                                                              "#B5DA77"), c("#DBBC74", "#B859DB", "#C792BF", "#96E37D", 
                                                                                            "#A0D6D1"), c("#BD55D9", "#B19CD7", "#D8CC6A", "#A8D7CF", 
                                                                                                          "#87E285", "#DC8085"), c("#BD55DA", "#B198D5", "#C1DD57", 
                                                                                                                                   "#DF6D83", "#9ED3DA", "#83E2A3", "#D8C093"), c("#B18DD2", 
                                                                                                                                                                                  "#7BE398", "#CFCF9A", "#7FD1D9", "#C8D851", "#BE51DB", 
                                                                                                                                                                                  "#D6C8CE", "#DB6F72"), c("#80D4DC", "#92E464", "#D6D0B7", 
                                                                                                                                                                                                           "#B947DE", "#DBC85C", "#8AE0A8", "#ACA6D5", "#C17BD0", 
                                                                                                                                                                                                           "#DD7377"), c("#CAD399", "#D888C6", "#7FDFC8", "#7F78D6", 
                                                                                                                                                                                                                         "#DACBC9", "#DBD14F", "#BE48DE", "#88B5D7", "#DE7665", 
                                                                                                                                                                                                                         "#83E474"), c("#DAD85A", "#D977CB", "#85C6D8", "#7D81D6", 
                                                                                                                                                                                                                                       "#D19A66", "#77E1BE", "#DC5F73", "#B947E1", "#D3DEBA", 
                                                                                                                                                                                                                                       "#86E271", "#D1AFCD"), c("#CAA7D0", "#D19965", "#DC6073", 
                                                                                                                                                                                                                                                                "#BA47E2", "#DBD650", "#D974C9", "#82E46F", "#7EC3D8", 
                                                                                                                                                                                                                                                                "#DADACC", "#C6DB99", "#7AE2C2", "#787ED6"), c("#8EE261", 
                                                                                                                                                                                                                                                                                                               "#CFDFAA", "#DBCECE", "#76E0A4", "#DCD554", "#D19962", 
                                                                                                                                                                                                                                                                                                               "#7FDCDA", "#DA5DD0", "#A344EA", "#7A74D3", "#DB6077", 
                                                                                                                                                                                                                                                                                                               "#D39AD3", "#7FA9CF"), c("#D39756", "#DA5DD0", "#7873D2", 
                                                                                                                                                                                                                                                                                                                                        "#CCE1D1", "#A344EA", "#DAD84F", "#D39AD4", "#D0D997", 
                                                                                                                                                                                                                                                                                                                                        "#84ADD2", "#DA5E75", "#D8B5B4", "#86E367", "#78DEA6", 
                                                                                                                                                                                                                                                                                                                                        "#71DDDF"), c("#DB47D3", "#7877D4", "#D2AED5", "#74E1A5", 
                                                                                                                                                                                                                                                                                                                                                      "#7CADD2", "#E0D250", "#D06692", "#D3DCD0", "#75DDDC", 
                                                                                                                                                                                                                                                                                                                                                      "#8FE25F", "#E4665A", "#DA85DB", "#CBDD9E", "#CA9B71", 
                                                                                                                                                                                                                                                                                                                                                      "#9545E8"), c("#9DE153", "#83DDB1", "#CDDB98", "#69E788", 
                                                                                                                                                                                                                                                                                                                                                                    "#D4E2D6", "#DB49D3", "#9645E8", "#DD5F73", "#72DCE0", 
                                                                                                                                                                                                                                                                                                                                                                    "#7777D5", "#D69A54", "#DED650", "#D1ABD7", "#7CACD2", 
                                                                                                                                                                                                                                                                                                                                                                    "#D77ECE", "#CCA69A"), c("#C3DCC5", "#DFD99D", "#AFDE88", 
                                                                                                                                                                                                                                                                                                                                                                                             "#E65566", "#D5A2E1", "#79D8E2", "#A144EA", "#CA729F", 
                                                                                                                                                                                                                                                                                                                                                                                             "#CA9781", "#DC5DD3", "#73E268", "#71DFB4", "#7BA3CF", 
                                                                                                                                                                                                                                                                                                                                                                                             "#7871D2", "#D99E49", "#DDCEDB", "#D8DC50"), c("#8861D8", 
                                                                                                                                                                                                                                                                                                                                                                                                                                            "#77E5D7", "#D1B0D9", "#9939EC", "#82C4DE", "#78E8A0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                            "#D7DB51", "#759E77", "#85E15B", "#DC47D3", "#DC88DD", 
                                                                                                                                                                                                                                                                                                                                                                                                                                            "#D06593", "#DE695B", "#738AD2", "#D5E3D7", "#D4E3A3", 
                                                                                                                                                                                                                                                                                                                                                                                                                                            "#D9A755", "#D6AC9F"), c("#A343EC", "#869BD8", "#BDE6CF", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#E09C48", "#E0576E", "#68E3D9", "#BBE68E", "#DBDA50", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#83CAE2", "#6DE4A2", "#79938B", "#7CE056", "#DCC9DC", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#7F68CD", "#D8978A", "#AEAC6B", "#D68ECC", "#E7E0B8", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#DE58D2"), c("#7DE2B2", "#DF7FD1", "#85C5DF", "#91E144", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#CCA1D3", "#9745EB", "#6AE1DD", "#AFE086", "#DC5C72", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#809474", "#DC47D2", "#E1D999", "#DED750", "#6C95D5", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#8168CB", "#D9964F", "#DED9E4", "#5AE888", "#D9A99E", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#CDE5CA"), c("#D3D8E5", "#76D9E3", "#DF7FD5", "#DEB2A9", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#E9556A", "#B96B8C", "#9645E9", "#7873D4", "#D7E74F", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#DEC44E", "#809874", "#D2A9DF", "#79DC4C", "#E3DC9E", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#CCE6CC", "#60EB94", "#77A4D0", "#DB47D2", "#D49057", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#78E3BC", "#B1E289"), c("#D1E199", "#DCDB50", "#DCAA4E", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#66E688", "#769B76", "#6C95D4", "#8B45E8", "#CDA9DC", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#B76F8D", "#D5A184", "#87C5DF", "#DFCED7", "#6BE1DD", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#D74DC5", "#D0E4CC", "#E56250", "#7E6BD0", "#97E254", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#DF89DC", "#83E6B6", "#EA5B8C", "#DE48EE"), c("#D9E3DB", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         "#E09852", "#DE89DB", "#A69A65", "#D6B1DD", "#7D6ACF", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         "#DEAFA6", "#6EA485", "#94B5D2", "#DED84F", "#8E47E8", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         "#E9576A", "#B0EBC0", "#E5E0A5", "#64EA9A", "#D84DC2", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         "#76D1E3", "#E147EE", "#B86A8B", "#ACDC7D", "#83E14D", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         "#67E8D7", "#7095DA"), c("#83CBE5", "#B4E78F", "#7CE758", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "#E8649E", "#7D95AD", "#BDE9CA", "#DAD8E3", "#6A90DD", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "#7BA27A", "#D69D4D", "#E4DAC1", "#DFDB4C", "#DE86DD", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "#E3DC96", "#6BE9A6", "#8565D2", "#6AE4DB", "#A56386", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "#963FEC", "#D1AADC", "#D79D8F", "#E8605A", "#85B741", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "#DB48D5"), c("#DFB8B1", "#DE9A3D", "#DDD8E9", "#DC46D2", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#558F44", "#859B8A", "#7397D3", "#E080D7", "#59DF52", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#9745EA", "#E3DB50", "#60EB95", "#81E2B4", "#ABE64C", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#E1E1A8", "#D18E6E", "#7C6DD0", "#E9556A", "#CABB67", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#83C7E0", "#BE6E94", "#B2E388", "#6AE3DF", "#CCE7D0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#D2A9DF"), c("#DB7EDB", "#E55760", "#DC46D3", "#E2D647", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#6E99DD", "#DAE2DE", "#9645E9", "#7EAF41", "#E4E1B2", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#D06996", "#D9C4E1", "#DCC675", "#74A17C", "#5BE986", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#BEEF66", "#B7E493", "#B9E7C9", "#69E1DD", "#76E8B5", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#7A6BCD", "#76E137", "#D98E4A", "#D5A4DE", "#8B839B", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#D7A89A", "#85C7E1"), c("#D76FD6", "#D2B6DC", "#7C9C9F", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#ACE083", "#D7A89F", "#E4615B", "#9631EB", "#E145DB", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#E2DC8B", "#C57496", "#B1E850", "#6993DC", "#7C5DAC", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#749D6B", "#D89A52", "#82C8E6", "#D5E5E5", "#61EB93", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#D69BE5", "#8862E5", "#E8569C", "#78E7B7", "#E0D6B3", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#6BE2DE", "#DDD346", "#57DD51", "#BFE9C3"), c("#D97ADC", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "#E3DD9B", "#80CAE1", "#905DE2", "#E8605D", "#E597CE", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "#78E04B", "#778698", "#D6E4E3", "#7CE2B2", "#D8EA51", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "#6283E6", "#DF45D7", "#D89652", "#825BA7", "#6AE3DE", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "#74A2DC", "#DC5793", "#C3A4E3", "#952FEB", "#C17986", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "#DAC3E1", "#DBC64B", "#859A71", "#ACDE81", "#E0B8A9", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "#C7E6C3", "#5DEB92"), c("#DF47EE", "#DFC574", "#DDDA4C", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "#D643BE", "#8BE33D", "#67E4DF", "#82CCE3", "#D9923D", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "#5AE784", "#CAE69E", "#E699BF", "#70A0DE", "#78E8B4", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "#7278D5", "#955ED8", "#C5A3E2", "#8539EA", "#DBE1E1", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "#E75E5B", "#8292A3", "#A45C7F", "#A6DF72", "#759F75", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "#E85691", "#D39D89", "#E2D9B6", "#DBC3E0", "#E083DA", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "#B8E7CA"), c("#996586", "#D7DC65", "#5BE57F", "#DF5993", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#D69C89", "#DF9A48", "#DD83DC", "#C0E4D7", "#E79ABF", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#BDECB8", "#76A1A4", "#CDA9E5", "#E4DCC2", "#A6DE7C", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#65A680", "#E0CFDD", "#83D1E9", "#D94BCB", "#E75D5C", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#E9E09C", "#96E542", "#7D6BD0", "#5D90D8", "#72E8B5", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#E3D43A", "#67E5E0", "#D337EF", "#A6A366", "#98B3DD", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#8B49E7"), c("#7377E8", "#E35F9B", "#E1BDDE", "#C7E6C6", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#DD99E2", "#D8E5E5", "#DA8F70", "#E95A5C", "#91E43B", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#E3D84A", "#6294DA", "#E241DA", "#9544EA", "#7D56A6", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#DF9C41", "#D66CD7", "#A9A170", "#ADE39C", "#BEEE7C", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#61AA84", "#5BE57F", "#6CEADF", "#E0BBAF", "#84C0E4", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#AD9FD5", "#B7708B", "#9BB44E", "#78D4E1", "#72E8B5", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#7C969D", "#EBE29E"), c("#E3E3CB", "#7DEAB0", "#716AD6", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#E9E19B", "#839494", "#62A982", "#7BBFE2", "#5BE783", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#CD6CD8", "#DCED53", "#A6A366", "#D498E4", "#E598C1", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#D98E70", "#8F40EA", "#B15B7E", "#C7B1E1", "#BDEABF", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#B2DADE", "#E3D7E9", "#8171A5", "#AADE7B", "#83E240", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#6898E1", "#5FE7D1", "#6FDDE2", "#E044DB", "#DE9942", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#DDB4AB", "#EB5AA4", "#D7C747", "#E8595D"), c("#83C3E4", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#7D50DA", "#E2E0BC", "#58E88E", "#D69BE5", "#D437B6", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#DCDFE9", "#E95169", "#E14EED", "#912DEC", "#A1AD6A", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#D5D960", "#5AC397", "#E677D0", "#DA876E", "#BCE6CF", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#63EDDC", "#A9E174", "#EDE397", "#7E9C98", "#7B96CD", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#E2B6B3", "#A5EBAE", "#667CD7", "#D9B1DB", "#BE6D94", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#D2ED41", "#76D9E1", "#E3C73B", "#B49976", "#76E348", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#A76CD8", "#DF9A46"), c("#98A880", "#83E639", "#599044", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "#5AE169", "#805AA5", "#7B8697", "#D6E5E6", "#DFC3BA", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "#902EEA", "#E5E6BB", "#DAE64D", "#E59B96", "#6BE9E6", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "#E7DD89", "#E292C3", "#ACE470", "#77A0DB", "#EA5769", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "#E045DD", "#83C9E3", "#64BBAD", "#DE69CA", "#5EEA97", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "#B6E7C8", "#D29DEB", "#B3E397", "#5C83E6", "#B35E80", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "#DCBF42", "#D2BAE0", "#935FE2", "#DC9152", "#B79A6E", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "#77E6B7"), c("#DFACA7", "#AFE687", "#78E8B8", "#E180DF", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#985CEB", "#E6A1DA", "#BE9E6B", "#912AEA", "#5CA175", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#57E356", "#5F6FD9", "#B6E8C3", "#A2E648", "#D0BCE0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#DD9346", "#E8635E", "#60EB93", "#81C8E5", "#66E1DE", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#E5D57D", "#A487D2", "#E5DEBD", "#B3DBD8", "#E75898", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#E246DC", "#DDC43F", "#929E86", "#B56889", "#6B99DD", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#7BAD3D", "#E8E9EB", "#E4EC55", "#A048B2", "#788A9E", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#D2E5A3"), c("#CC7293", "#EAE19A", "#E95A5C", "#BCE9C0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "#DF78E0", "#7E9591", "#985BE4", "#A48CAF", "#925E56", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "#E4E3CA", "#5BE680", "#ACE081", "#8431EA", "#A7BAEA", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "#5BAA84", "#D33DBB", "#DDAFA4", "#E95595", "#6DD4DC", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "#5D92D8", "#D9DC4F", "#71BFDF", "#76E8B2", "#E18DD4", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "#E3AF3C", "#A7A567", "#5974E5", "#E146EF", "#7D56A4", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "#DF915E", "#B996E7", "#69EDE1", "#92E442", "#E2B2E1", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "#B4DCE1", "#E6D7E5"), c("#60B94B", "#7BE533", "#97A46F", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#E3DFBE", "#A05A7D", "#962DEB", "#5EE8D3", "#DF44D7", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#E793C5", "#71DDE2", "#E4AD3A", "#B7E6C6", "#80EBB3", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#DBD565", "#E95A5C", "#D3B172", "#E5EBA9", "#AFD6DF", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#6595DA", "#9659E3", "#746ED0", "#7C9491", "#DBC4E6", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#E65393", "#D79489", "#D1A4EC", "#B2E766", "#A593B9", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#EAEFED", "#B1E390", "#DDBBB5", "#D88752", "#5AED8C", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#DFE13D", "#D675D6", "#7ABFE3", "#54AC85"), c("#E2CFBD", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#B9A2DF", "#9559E3", "#748CA2", "#814E6C", "#902DEB", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#DD7FDF", "#DFD944", "#C1F174", "#E6DC88", "#DFA342", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#B1E394", "#9ECCE9", "#6597DD", "#D9C2DF", "#94EDB8", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#CF739F", "#61C1DD", "#EBA6E0", "#81B53F", "#A5A167", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#5BEACC", "#DEE9EA", "#70E2E4", "#DD8D70", "#DAA8A3", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#5573E5", "#58AE83", "#E4543D", "#5BEA8A", "#E4E7B9", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#8FA99B", "#E44DEB", "#EA547C", "#BBE7CD", "#D33DB4", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#83E539", "#8A69BF"), c("#D643BF", "#9AB9E4", "#D89078", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         "#B6D8D3", "#7C918C", "#E1DCBA", "#DFC5E5", "#A38CB0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         "#5E91D7", "#89E346", "#9375CC", "#A7E07E", "#663797", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         "#5677E8", "#E7BD3E", "#69EBE7", "#65C1B4", "#EAECEC", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         "#8133E9", "#D8E444", "#CDA7EA", "#E95A5C", "#E75592", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         "#D58A3F", "#DFB7B0", "#EB9AC4", "#73E8B5", "#B7EBC4", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         "#A561E9", "#A85C7E", "#DCE9A4", "#73A068", "#E045EF", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         "#99DAEF", "#C8AF72", "#58BBDA", "#5AE784", "#DBD769", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         "#E183DB"), c("#E0B9B2", "#E45B97", "#B4DB96", "#6181E7", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#5FE8CD", "#E9E754", "#D9C3E4", "#71A1DD", "#975EE0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#EBE592", "#A4CEE3", "#5AEC94", "#94EDB6", "#DD9679", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#64BED6", "#B0EA7F", "#E64FEC", "#EA595B", "#6FE34A", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#86849D", "#6833E6", "#E2E9E8", "#71E5E8", "#E095BF", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#B46176", "#DE84E0", "#D33EB8", "#B333EF", "#C1A970", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#E2E1BB", "#53A53D", "#8059A5", "#B7E4CA", "#82957A", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#EBB93F", "#57B189", "#BCEC4C", "#D3883D", "#C9A3E5", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#B7BB44"), c("#BDA46F", "#987EA4", "#6294DB", "#D88A69", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#DA7EDE", "#59AE45", "#ACB9BD", "#E046EF", "#E0ACA6", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#9EBAE6", "#BEE8C5", "#E65593", "#E3EF54", "#E8D5E4", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#5CEC9A", "#56C199", "#5FEA71", "#D1EBEB", "#80D1EA", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#99AD8D", "#5E7F4A", "#E95A5C", "#84F0BF", "#985BE4", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#5F9EA9", "#81E537", "#E8C647", "#6AE6E1", "#726ACD", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#D644C0", "#E0B6DF", "#EAE394", "#8431EA", "#B5BF43", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#E790C5", "#E7DEBF", "#AD5C7D", "#C8A2EA", "#D99138", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#B5E296", "#AEEB76"), c("#E241D9", "#64EDB8", "#786BD0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#7DBCE5", "#E1CFAE", "#E3E8E4", "#ECDF8E", "#E1E834", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#E95A5C", "#5ED2E1", "#9E5B7C", "#E8C544", "#9443EB", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#B19565", "#AA9ED1", "#A7D8EA", "#D06DD8", "#E497BA", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#6EE330", "#B5E9C7", "#B6BD42", "#D69033", "#E45A9D", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#61834A", "#5EBF95", "#ACEBA9", "#B6EC84", "#65BE4D", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#D999E3", "#62EADF", "#9EB39E", "#E39370", "#E9F272", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#E6EDBE", "#9FE8E2", "#5DED8A", "#DDC4E4", "#ACC27D", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#6391D8", "#DEB1AC", "#6E9DA5", "#A8E751"), c("#EDB83F", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#EBDE91", "#9BB8E4", "#6FAC3C", "#D33DB9", "#6CC5E4", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#DEC3E3", "#AB6486", "#E5E73B", "#AAA462", "#A4EB7F", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#DF4AED", "#E997D0", "#CBC547", "#E14B3E", "#8B43E9", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#B2E099", "#55E455", "#E56D7C", "#B8E7C8", "#6293DA", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#7175E8", "#D8A7AB", "#E75595", "#A2E443", "#D87CE0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#EAEEEB", "#97997B", "#E2CAB8", "#C3A2E1", "#778A9C", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#59EB96", "#67BFBA", "#5CE8D3", "#D79048", "#DC9579", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#7D55A6", "#78E9ED", "#53A87A", "#85EBB7", "#E4E7BB", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#B1D7E0", "#DDF17C"), c("#CF32AE", "#E0ED55", "#E64B44", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "#80B23C", "#DAB473", "#5BE986", "#E7E9EB", "#6164E3", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "#E9547D", "#A1D4E9", "#757699", "#5B89E4", "#BAE3D2", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "#9E9C73", "#6B27E7", "#7BA9E0", "#5BC1E0", "#DCC646", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "#D78F35", "#E479DB", "#79A2A3", "#7A56A1", "#E2E394", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "#E25EE7", "#B1E783", "#6CE4E6", "#DDB4AD", "#E837EA", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "#DAC3E2", "#E4DFBC", "#E16EAB", "#9E57E2", "#DA856F", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "#C0A6E3", "#85E541", "#B439EF", "#E198E1", "#80E8AB", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "#E29ABE", "#62E6CA", "#BFEDBB", "#AB5F7C", "#64A87A", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "#B07CE2"), c("#5BE37C", "#EA595A", "#DAD664", "#E45996", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#7F4E9C", "#DEC5E4", "#D63BB6", "#DCA237", "#739BA8", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#61C4C3", "#E2DEBE", "#777CC9", "#E1DE3B", "#74EDF0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#579866", "#EB99C3", "#386FDF", "#B0EEBE", "#C1EAD8", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#9A7CEC", "#87D1EA", "#E03BEC", "#842CE9", "#EAEFAD", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#B86376", "#ADB46F", "#96B7E5", "#AEE264", "#63E5AB", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#AD8F6E", "#8E54E6", "#CCA6EA", "#65EAD1", "#ACE491", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#AF8DB2", "#E6C883", "#E28E62", "#D164E5", "#6B6F90", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#5999E0", "#DEE5E9", "#82E639", "#9BB89F", "#E185DC", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#E3B6AF"), c("#69B6D8", "#5EED9C", "#E035EB", "#DB7EDF", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "#B0E496", "#8654DF", "#A3EEBF", "#DFB3AE", "#C1E6D3", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "#82D9F0", "#E2B0D7", "#B2D3E0", "#D2C543", "#E579B6", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "#5C8D3D", "#DFEB45", "#DB9D37", "#5AEACB", "#E3D8C0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "#D73AB5", "#81E537", "#5ADF69", "#757DEB", "#6EECEB", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "#D4BB74", "#EB5080", "#6096DD", "#A6B580", "#D15FE7", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "#E27F82", "#7F7A9E", "#B0688C", "#AEB7E8", "#852EEA", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "#7D9795", "#E6EAB8", "#7F5BA9", "#EDED8E", "#5FBB90", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "#65C4C1", "#D59FE7", "#B0E66E", "#E54C41", "#E9E0EC", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "#E49161", "#AB8668"), c("#97DE94", "#5CEF79", "#6994A7", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#60C3C2", "#B19170", "#E8F071", "#EBCA67", "#C2A1E9", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#59BE95", "#EA595A", "#5CECD9", "#E1C1E3", "#D43FBC", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#84E637", "#57EB9A", "#6C74E7", "#DCE1C1", "#EEE7A9", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#7D55A5", "#C6E591", "#8B43E9", "#D87BE2", "#C16680", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#BBEDC0", "#75E9F0", "#ECEEEF", "#E2DD34", "#D5912E", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#B7E8D9", "#AFED77", "#E4AAA5", "#6393DB", "#ABBE3F", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#E95595", "#E39268", "#82F0BE", "#5ABA49", "#E492CD", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#E2C9BD", "#6DC8E8", "#9FB9E6", "#DF49EE", "#B6B370", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#94AC95", "#AFD5E2", "#A185A5", "#5F8449"), c("#D333B3", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#65B796", "#97819D", "#4A81B9", "#E8548F", "#D5B972", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#BDE095", "#E89BC0", "#D55DE6", "#EA5859", "#7DE439", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#E3DEBE", "#8DBAE5", "#5FEDE1", "#A7A17C", "#5BE984", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#78529F", "#709AEB", "#9AE38E", "#9AEFBE", "#66C7E4", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#E678DB", "#77DEE2", "#D18738", "#BDE7C9", "#4E71E3", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#B3DAE3", "#EAECEB", "#E2DD65", "#AE7BE4", "#B2A0DB", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#ECEDAA", "#5B814B", "#B76175", "#8028E7", "#DDC4E5", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#DE9273", "#E035EB", "#7AAC3B", "#E1B7B0", "#B6ED6F", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#E2A5EE", "#D177B8", "#EAB83C", "#789EA2", "#52E4AF", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "#9155E8", "#D8DF32"), c("#6CEBE4", "#E7A09E", "#5E83DE", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#E9609A", "#BF9EEA", "#EAEFE7", "#B0DF99", "#E9B2DA", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#DDD4ED", "#65CFE2", "#6391A5", "#D53EBB", "#D48F30", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#F1CE78", "#A7E458", "#E33DED", "#DEBEB9", "#852FE9", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#8155A7", "#E95A5C", "#E89CE3", "#7467E7", "#B19170", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#90E331", "#59C4A0", "#B5DCD5", "#B65CE9", "#AC9DC5", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#B1BD44", "#A35869", "#B7709E", "#E3DDB9", "#5AE169", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#58E530", "#578E42", "#BCB575", "#86A58E", "#83F0BF", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#B2EB84", "#EAEFAB", "#5CEA97", "#9DCEE6", "#E2EC35", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#E39268", "#DD7ADD", "#6DA3DE", "#E3C63A", "#B9EAC6", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "#E7F074"), c("#C3A4E8", "#5693D2", "#E3DD3D", "#DED662", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#E75CE9", "#B4CFD1", "#DC9939", "#DFD1B9", "#99D7EE", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#5BE66E", "#9AE9E8", "#B2DD97", "#61E6C4", "#E6EBBC", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#E6B2AE", "#5AA87B", "#5BC1C0", "#7783D2", "#59BDDE", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#B863DD", "#B5EB54", "#EA5491", "#9AB8E5", "#B7E5C6", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#5DECE9", "#D12FAF", "#D88766", "#EDDA91", "#7E4FDD", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#C76E8D", "#7A9592", "#80509C", "#93EDB4", "#E4B3DD", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#862FEC", "#967BEC", "#3E7CE9", "#E377D7", "#A8A969", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#72E22D", "#E1D0E6", "#76AB3A", "#DA37ED", "#EB5758", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#9F86A6", "#E291DA", "#B29474", "#EAF1EC", "#58EB98", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "#B3EB83"), c("#A8E6D9", "#E5503D", "#D6E2C5", "#DCD1EA", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "#AEBF7F", "#B29270", "#76E22F", "#63ECE0", "#EDECA5", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "#DF8665", "#5DE870", "#E09BE4", "#BBEDC2", "#E35EA7", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "#85EEBB", "#985AE6", "#E85376", "#618C9F", "#E2E533", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "#5072E3", "#EAF0EE", "#8F2CEA", "#A3E783", "#57EC9A", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "#E6B9E3", "#82A590", "#DAB773", "#EDE36D", "#7EBCE5", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "#E044DD", "#E4C5BD", "#B5EB61", "#764D9B", "#5ABF95", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "#E8A29E", "#AC87E2", "#A6627E", "#B0CBCE", "#B59DA8", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "#99DAEE", "#6393DA", "#B7BD3C", "#D96ED5", "#E690BB", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "#61824D", "#56A83D", "#DCA138", "#B0A1D6", "#BBE898", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "#5BD0DD", "#EBDCB5"), c("#AFEB76", "#7BA29F", "#E08B64", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "#5CEB92", "#ADE3C8", "#E33CEC", "#5EAF84", "#DCF1B3", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "#EDD643", "#EDEC89", "#ACE091", "#D2902F", "#AE9372", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "#4F81EB", "#8B5651", "#CEE6E8", "#BAA8E8", "#6DE1E4", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "#EADDA9", "#B3B374", "#E4D9C1", "#E29DE7", "#7A6AB6", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "#E7512D", "#59E7D2", "#E8A7A0", "#B177A7", "#988EAE", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "#7ACDE7", "#D7ECCC", "#E3B6DD", "#BEBD42", "#BD5EE9", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "#852FE9", "#90BAE5", "#793594", "#CFED49", "#69993C", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "#B07CE2", "#BE6183", "#E95F68", "#DABCBB", "#F08EBC", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "#D83CBD", "#6CDF48", "#E85091", "#5D93D5", "#735EE6", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "#F1BD76", "#E278D8", "#E0D5EC", "#89EEBA"), c("#7BB039", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#81E537", "#DE6CAA", "#67ECE1", "#E45730", "#E84986", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#CE32AF", "#6AE7B7", "#E2B4B0", "#5DEB97", "#F0BD77", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#E295DF", "#B2E394", "#EBE59D", "#8565AF", "#E0C338", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#DFEB4B", "#E1D6C7", "#9CBAE7", "#D39031", "#DDCEA9", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#AF7BE6", "#B1D9E6", "#5178E9", "#B1ED79", "#4DA67B", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#A28CB0", "#6BD9E2", "#ECEFF1", "#BEEADA", "#6296A9", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#6194DB", "#7F28E7", "#C9A8EB", "#DFC5E5", "#AE6279", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#8DC0B4", "#E85F69", "#A7BC89", "#DE5DE6", "#A39859", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#E1EEC4", "#D98D6F", "#EC9EC3", "#62369D", "#5BE66E", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#E479DB", "#6FC7E7", "#9758ED", "#DDDB6B", "#AAEEBF", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "#918C7D", "#DD34EB"), c("#DFA031", "#ECEEEF", "#E043D5", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#E9DBAA", "#9273CA", "#EA4F83", "#B1A3D9", "#92EFC0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#9B5BE8", "#B0D9E8", "#DC8753", "#6FE22D", "#70DCEE", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#F195B5", "#5AEB96", "#B0587D", "#B4DDD2", "#6BEAE1", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#A27898", "#53B6BD", "#DC9EEA", "#7D3693", "#5F92D7", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#A7AD8F", "#A7E652", "#697E48", "#B5EE84", "#952CEB", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#EBB1DF", "#E0E7C9", "#7B989F", "#58B388", "#E64942", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#81BFE6", "#EBE388", "#A0E095", "#5BE66E", "#E2ADA6", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#5273E5", "#E7F470", "#DF7877", "#E5CCC0", "#DBC6E2", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#B6E8C4", "#E177BC", "#7AAC3B", "#D9C843", "#AC8567", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#E1EB34", "#B2C27F", "#DBF2B2", "#DD76DF", "#D9B06C", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "#50E7B8"), c("#558F44", "#998DB0", "#B27DA5", "#E295E1", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#E3D2E9", "#8E2CE9", "#76E7EC", "#E4453F", "#D76AA4", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#4E71E3", "#9777D0", "#E8EBB3", "#E5879E", "#DC79DF", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#E54DED", "#EBBA3E", "#E3D9BD", "#EDF0EE", "#713E8F", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#66C7E5", "#81B9E4", "#7CE431", "#72DFAD", "#EDE987", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#D7BBBA", "#ABBC7A", "#ACBAE6", "#DDB9E3", "#E5EA43", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#9A5DE9", "#5CBBB1", "#B2E96A", "#4A82BB", "#6C8E94", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#C936B6", "#D7B974", "#C0EAD8", "#D28946", "#AED5E2", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#B6EEC0", "#975564", "#AB8F6E", "#EB507F", "#EEA5D2", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#E954B9", "#6E9CED", "#C5A5EB", "#BCC03F", "#5AE068", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#E87F73", "#5EEDDA", "#9BBAA0", "#5AEB97", "#E9ADA0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "#ACE591"), c("#B3EC80", "#E4AD3A", "#7B54A0", "#AED48C", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "#4880B8", "#83D9F0", "#EE76AE", "#AC6FE5", "#DD30ED", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "#C2EEBF", "#8129EA", "#64C399", "#ECECEC", "#6EB9E0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "#C36BAB", "#C1A7A8", "#DAE444", "#5FECCD", "#A059F0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "#CEADE5", "#A0B8E8", "#E99FE8", "#EA5678", "#DF3CA0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "#B86673", "#E69CC0", "#58EB98", "#5AAAB4", "#518B58", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "#79AC3B", "#BDE4D1", "#DDD764", "#E2513C", "#E02CBE", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "#6D9BEB", "#4FE1E4", "#87E63D", "#8AE8E1", "#94F0B5", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "#5C40C4", "#DE9262", "#E1DCBA", "#EBB7AB", "#5BE66E", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "#EF6FD4", "#DEC6E3", "#917FA1", "#EDE19A", "#AF8EE0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "#B2D3DF", "#5277E8", "#E085E8", "#B29F66", "#859685", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "#E15DED", "#BA41B5"), c("#E0A49A", "#627842", "#BBE7C4", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#99EFBE", "#8BE9E4", "#5FE6D0", "#E95493", "#EFEDEF", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#7CAD3D", "#5AEC94", "#B49C64", "#B6D4DE", "#E3E3C4", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#8E5B53", "#57BCBE", "#E49268", "#DFC4E5", "#DDC0BD", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#C9A4E9", "#B07AE7", "#A389AD", "#E033EB", "#52BADD", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#88D9ED", "#98BFE2", "#862DEB", "#B2C580", "#EAEEAD", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#ADEB6D", "#73919B", "#C56B8F", "#4E71E3", "#E95A5C", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#E1CCA0", "#58B184", "#5BE66D", "#95AA91", "#ECCF75", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#7E59A4", "#D48F30", "#4780BA", "#94B2E9", "#DAC63B", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#E7EE70", "#9BE38D", "#D55EE9", "#56E7B1", "#76E22F", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#C5EFA1", "#8650DF", "#D439B6", "#E27EDB", "#C8EADD", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "#EA99D1", "#4DEAF0", "#6E99EB", "#D5EC37"), c("#B29171", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#AABBBF", "#BFA4E7", "#92A88D", "#AEDAEB", "#98F0BC", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#D13BBA", "#DDF2AF", "#5DEC94", "#BFE4D8", "#CB70B0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#E3CFBA", "#B0C17F", "#E7567E", "#E7DAA2", "#4EA6D3", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#E28870", "#BDF274", "#E0CFE7", "#80E435", "#9677D1", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#EDF0EF", "#E247EE", "#4F72E4", "#698996", "#E64B44", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#6F418F", "#5BECD9", "#A78CAF", "#4FAE7E", "#E2E8C7", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#6393DC", "#DBC540", "#B9E9C5", "#9BB9E6", "#DAB069", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#892DE9", "#DA8E33", "#A95F74", "#4FE9F1", "#EC95B9", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#78C8AE", "#81B53F", "#E3B8E2", "#56BCC3", "#995DE8", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#DD78E2", "#E49CE7", "#5BE66D", "#677B44", "#8BE8E7", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#54E5B2", "#E1E942", "#A8E38E", "#F063B6", "#E3B3AF", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "#E9E57E", "#75D0EE"), c("#E955EA", "#58A33C", "#AED3DD", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#6116BC", "#D8F2B2", "#5BE66E", "#A9E58A", "#B4BF3F", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#E5AFF1", "#A24678", "#DFCAE8", "#EDE988", "#9E9777", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#D59031", "#EB656D", "#E5EA43", "#DB2BAF", "#61EBDB", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#DF4432", "#7611F5", "#6E8B99", "#EB5391", "#9FE9E5", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#8CBAE6", "#B2EBC7", "#5178E9", "#8263AF", "#E17E98", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#76E22F", "#E7CCC3", "#ED75DB", "#7AE6B2", "#AF6EE7", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#E7C542", "#5F93D9", "#57EC9A", "#61AC8A", "#AFA0D2", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#DC8BE6", "#B6EC65", "#7136F5", "#8D57ED", "#E8AD9D", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#59D7DE", "#B096EA", "#DC82C3", "#9C667D", "#D88665", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#C935F0", "#BB41B8", "#5F39A2", "#D5E2C6", "#E6DAAA", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#D9B36C", "#E4A5CC", "#C0A3AA", "#A5B875", "#EBF0EE", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "#68C8E6"), c("#51B5D8", "#6AEFE7", "#83E639", "#D28064", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#66E5BF", "#639ADF", "#D4E3F0", "#E66E81", "#52E468", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#89569F", "#EDA5D2", "#7376C0", "#8C56ED", "#EAF1E7", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#9B7BEB", "#B7EDC8", "#C2F2AB", "#70EA89", "#417DE9", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#E4CEBC", "#C632EF", "#98BEEA", "#B19474", "#ECADA2", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#5BC6C3", "#ED4CE4", "#E198E8", "#3DEC9F", "#AFA057", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#E0C8E7", "#EACB85", "#D9E835", "#6B27E7", "#C7AAB1", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#AE6A89", "#7DE0F1", "#5B38A2", "#D970B2", "#B0D8CB", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#D88B3F", "#C7A6EB", "#CD31AD", "#EBEFA4", "#E74886", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#A6CA83", "#93C5D6", "#C360E3", "#E67ADE", "#EDB83D", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#58B184", "#AD9ECB", "#DFE3B9", "#7AB642", "#8EEEB0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#7F9673", "#CAC442", "#B3EE7D", "#ECF074", "#6E7E8E", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "#E64B44"))
  l[[k]]
}


normalize.nQuantiles <- function (x, probs = 0.5, shareFeature = FALSE, ref = 1) 
{
  if (is.data.frame(x)) 
    x <- apply(x, 2, as.numeric)
  if (shareFeature) {
    if (length(probs) == 1) {
      fac <- vapply(seq_len(ncol(x)), function(i) {
        ir <- !is.na(x[, ref]) & !is.na(x[, i])
        quantile(x[ir, ref], probs = probs) - quantile(x[ir, 
                                                         i], probs = probs)
      }, numeric(1))
      x <- sweep(x, 2, fac, "+")
    }
    else {
      for (i in seq_len(ncol(x))) {
        ir <- !is.na(x[, ref]) & !is.na(x[, i])
        refquant <- quantile(x[ir, ref], probs = probs, 
                             na.rm = TRUE)
        colquant <- quantile(x[ir, i], probs = probs, 
                             na.rm = TRUE)
        mod <- lm(refquant ~ colquant)
        x[, i] <- mod$coefficients[[2]] * x[, i] + mod$coefficients[[1]]
      }
    }
  }
  else {
    if (length(probs) == 1) {
      colquant <- colQuantiles(x, probs = probs, na.rm = TRUE)
      grandquant <- quantile(x, probs = probs, na.rm = TRUE)
      off <- colquant - grandquant
      x <- sweep(x, 2, off, "-")
    }
    else {
      colquant <- t(colQuantiles(x, probs = probs, na.rm = TRUE))
      grandquant <- quantile(x, probs = probs, na.rm = TRUE)
      for (i in seq_len(ncol(x))) {
        mod <- lm(grandquant ~ colquant[, i])
        x[, i] <- mod$coefficients[[2]] * x[, i] + mod$coefficients[[1]]
      }
    }
  }
  x
}


normalize.totsum <- function (x) 
{
  x <- 10^x
  xs <- colSums(x, na.rm = TRUE)
  x <- sweep(x, 2, xs, "/")
  x <- x * median(xs)
  log10(x)
}


normalizeColWise <- function (x, method = c("Median centering", "Median centering (shared ID)", 
                                            "Total sum", "median centering + variance stablization")[1]) 
{
  nr <- which.max(colSums(!is.na(x)))
  normfun <- switch(method, None = function(x) x, `Median centering` = function(x) normalize.nQuantiles(x, 
                                                                                                        probs = 0.5, shareFeature = FALSE, ref = nr), `Median centering (shared ID)` = function(x) normalize.nQuantiles(x, 
                                                                                                                                                                                                                        probs = 0.5, shareFeature = TRUE, ref = nr), `Total sum` = normalize.totsum, 
                    `median centering + variance stablization` = function(x) normalize.nQuantiles(x, 
                                                                                                  probs = seq(0.25, 0.75, by = 0.05), shareFeature = FALSE, 
                                                                                                  ref = nr))
  normfun(x)
}


normalizeData <- function (x, colWise = c("None", "Median centering", "Median centering (shared ID)", 
                                          "Total sum", "median centering + variance stablization")[1], 
                           rowWise = c("None", "Reference", "Batch mean", "Batch reference")[1], 
                           ref = NULL, batch = NULL) 
{
  d <- normalizeColWise(x, method = colWise)
  if (rowWise == "Reference") {
    if (is.null(ref)) 
      stop("Reference not given!")
    ina <- is.na(d)
    d[ina] <- min(d, na.rm = TRUE)
    d <- removeVarQC(d, ref)
    d[ina] <- NA
  }
  else if (grepl("Batch", rowWise)) {
    if (is.null(batch)) 
      stop("Batch not given!")
    if (is.null(ref) && rowWise == "Batch reference") 
      stop("Reference not given!")
    d <- rowshift(d, batch = batch, ref = ref, useMean = rowWise == 
                    "Batch mean")
  }
  d
}


null2empty  <- function (x) 
{
  if (is.null(x)) 
    return("")
  x
}


# omicsViewer <-  function (dir, additionalTabs = NULL, filePattern = ".(RDS|DB|SQLITE|SQLITE3)$", 
#             ESVObj = NULL, esetLoader = readESVObj, exprsGetter = getExprs, 
#             pDataGetter = getPData, fDataGetter = getFData, defaultAxisGetter = getAx, 
#             appName = "omicsViewer", appVersion = "1.1") 
#   {
#     app <- list(ui = fluidPage(app_ui("app")), server = function(input, 
#                                                                  output, session, aTabs = additionalTabs, f_eset = esetLoader, 
#                                                                  f_exprs = exprsGetter, f_pd = pDataGetter, f_fd = fDataGetter, 
#                                                                  axg = defaultAxisGetter) {
#       callModule(app_module, id = "app", .dir = reactive(dir), 
#                  additionalTabs = aTabs, filePattern = filePattern, 
#                  esetLoader = f_eset, exprsGetter = f_exprs, pDataGetter = f_pd, 
#                  fDataGetter = f_fd, defaultAxisGetter = axg, appName = appName, 
#                  appVersion = appVersion, ESVObj = reactive(ESVObj))
#     })
#     runApp(app)
#   }


parseDatTerm <- function (file, outputDir = NULL, ...) 
{
  message("Reading dat file ...")
  d0 <- readLines(file, ...)
  d0 <- split(d0, cumsum(d0 == "//"))
  org <- trimws(sub("OS", "", grep("^OS", d0[[1]], value = TRUE)))
  org <- make.names(paste(org, collapse = ""))
  while (grepl("\\.\\.", org)) org <- gsub("\\.\\.", ".", org)
  fn <- basename(file)
  vn <- paste0("_", gsub("-", "", Sys.Date()), "_", org, ".annot", 
               sep = "")
  if (is.null(outputDir)) 
    outputDir <- dirname(file)
  outputFile <- file.path(outputDir, sub("(.dat|.dat.gz)$", 
                                         vn, fn))
  message("Processing ...")
  dd <- lapply(d0, function(x) {
    dr <- grep("^DR", x, value = TRUE)
    an <- stringr::str_split_fixed(trimws(sub("^DR", "", 
                                              dr)), ";", 4)
    if (nrow(an) == 0) 
      return(NULL)
    ac <- stringr::str_split_fixed(trimws(sub("^AC", "", 
                                              grep("^AC", x, value = TRUE))), ";", n = 2)[1]
    name <- stringr::str_split_fixed(trimws(sub("^ID", "", 
                                                grep("^ID", x, value = TRUE))), " ", 2)[1]
    gn <- trimws(sub("^GN", "", grep("^GN", x, value = TRUE)))
    gn <- grep("Name=", gn, value = TRUE)
    if (length(gn) == 0) 
      gn <- NA
    else {
      gn <- gsub("Name=|;$", "", strsplit(gn[1], " ")[[1]][1])
    }
    data.frame(ID = name, ACC = ac, geneName = gn, source = an[, 
                                                               1], term = an[, 2], desc = an[, 3], stringsAsFactors = FALSE)
  })
  dd <- do.call(rbind, dd)
  uid <- paste(dd$source, dd$term)
  tb <- table(uid)
  i <- which(uid %in% names(tb[tb >= 5]) & uid %in% names(tb[tb < 
                                                               0.25 * length(d0)]))
  dd <- dd[i, ]
  for (i in colnames(dd)) dd[[i]][is.na(dd[[i]])] <- "_NA_"
  if (!is.null(dd)) {
    message("Writing table ...")
    write.table(dd, file = outputFile, col.names = TRUE, 
                row.names = FALSE, quote = FALSE, sep = "\t")
  }
  invisible(dd)
}


phenoTemplate <- function (label, quant = c("LF", "TMT", "undefined")[1]) 
{
  n <- min(str_count(label, "_")) + 1
  tab <- data.frame(Label = label, stringsAsFactors = FALSE)
  if (quant %in% c("TMT", "undefined")) {
    cn <- str_split_fixed(label, "\\.", 2)
    tab$Batch <- cn[, 2]
    tab$Channel <- str_extract(cn[, 1], "(\\d)+")
  }
  else {
    tab$Batch <- "1"
  }
  tab$Reference <- FALSE
  if (n > 1) {
    m <- str_split_fixed(label, pattern = "_", n = n)
    colnames(m) <- paste0("Var", seq_len(ncol(m)))
    tab <- cbind(tab, m)
  }
  it <- which(vapply(tab, is.factor, logical(1)))
  if (length(it) > 0) 
    tab[it] <- lapply(tab[it], as.character)
  tab
}


plot_roc_pr_module <- function (input, output, session, reactive_param, reactive_checkpoint = reactive(TRUE)) 
{
  output$roc_pr <- renderPlot({
    req(reactive_checkpoint())
    draw_roc_pr(reactive_param()$y, reactive_param()$x)
  })
}


plot_roc_pr_ui <- function (id) 
{
  ns <- NS(id)
  shinycssloaders::withSpinner(plotOutput(ns("roc_pr")))
}


plotly_barplot <-
function (x, names, highlight = NULL, highlight_color = "red", 
          highlight_width = 1, highlight_legend = "highlighted", background = NULL, 
          background_color = "gray", background_width = 1, background_legend = "background", 
          ylab = "ylab", xlab = "xlab", sort = c("none", "increasing", 
                                                 "decreasing")[1], tooltips = NULL, source = "plotlybarchart") 
{
  background <- na.omit(background)
  highlight <- na.omit(highlight)
  data <- data.frame(x = x, names = names, stringsAsFactors = FALSE)
  sort <- match.arg(sort, choices = c("none", "increasing", 
                                      "decreasing"))
  if (sort == "increasing") {
    data$xpos <- rank(data$x, ties.method = "random")
  }
  else if (sort == "decreasing") {
    rk <- rank(data$x, ties.method = "random")
    data$xpos <- (max(rk) + 1) - rk
  }
  else data$xpos <- seq_len(nrow(data))
  fig <- plot_ly(data, source = source)
  fig <- add_lines(fig, x = data$xpos, y = data$x, type = "lines", 
                   line = list(color = "lightgray"), name = "Ranking stats")
  st <- (max(data$x, na.rm = TRUE) - min(data$x, na.rm = TRUE))/20
  if (!is.null(tooltips)) 
    txt <- tooltips
  else txt <- names
  if (!is.null(highlight)) {
    y <- data$x[highlight]
    y <- sign(y) * pmax(abs(y), st)
    fig <- add_bars(fig, x = data$xpos[highlight], y = y, 
                    type = "bar", name = highlight_legend, marker = list(color = highlight_color), 
                    opacity = 0.8, width = highlight_width, text = txt[highlight], 
                    hoverinfo = "text")
    fig <- plotly::layout(fig, xaxis = list(showticklabels = TRUE, 
                                            tickvals = data$xpos[highlight], ticktext = data$names[highlight]))
  }
  if (!is.null(background)) {
    y <- data$x[background]
    y <- sign(y) * pmax(abs(y), st)
    fig <- add_bars(fig, x = data$xpos[background], y = y, 
                    type = "bar", name = background_legend, marker = list(color = background_color), 
                    opacity = 0.8, width = background_width, text = txt[background], 
                    hoverinfo = "text")
  }
  fig <- plotly::layout(fig, yaxis = list(title = ylab), xaxis = list(title = xlab), 
                        legend = list(orientation = "h", xanchor = "center", 
                                      x = 0.5, yancher = "center", y = 1.1))
  fig <- plotly::config(fig, toImageButtonOptions = list(format = "svg", 
                                                         filename = "omicsViewerPlot", width = 700, height = 700))
  toWebGL(fig)
}


plotly_barplot_module <- function (input, output, session, ...) 
{
  ns <- session$ns
  output$plot <- renderPlotly(do.call(plotly_barplot, args = list(..., 
                                                                  source = ns("barplotly"))))
  reactive(event_data("plotly_click", source = ns("barplotly")))
}

plotly_barplot_ui <- function (id) 
{
  ns <- NS(id)
  shinycssloaders::withSpinner(plotlyOutput(ns("plot")), type = 8, 
                               color = "green")
}


plotly_boxplot <- function (x, i = NULL, highlight = NULL, ylab = "ylab", extvar = NULL, 
                            ylab.extvar = "ylab.extvar") 
{
  o.name <- colnames(x)
  colnames(x) <- paste0("B", str_pad(seq_len(ncol(x)), pad = "0", 
                                     nchar(ncol(x))))
  mat_i <- x[i, , drop = FALSE]
  tlpmap <- structure(o.name, names = colnames(x))
  if (length(highlight) == 0) 
    highlight <- NULL
  if (length(i) == 0) 
    i <- NULL
  convertDF <- function(m, c, maxr = 200) {
    if (nrow(m) > maxr) 
      m <- apply(m, 2, quantile, probs = seq(0, 1, by = 0.02), 
                 na.rm = TRUE)
    df <- reshape2::melt(m)
    hp <- colnames(m)[c]
    df$Cat <- c("n", "c")[as.integer(df$Var2 %in% hp) + 1]
    na.omit(df)
  }
  df <- convertDF(x, c = highlight)
  fig <- plot_ly(showlegend = FALSE)
  fig <- add_boxplot(fig, data = df, x = ~Var2, y = ~value, 
                     color = ~Cat, colors = c("#35c636", "#8a8588"), boxpoints = is.null(i), 
                     opacity = ifelse(is.null(i), 1, 0.1), text = paste0("<b>Sample: </b>", 
                                                                         tlpmap[df$Var2]), hoverinfo = "text")
  if (nrow(mat_i) > 20) {
    dfi <- convertDF(mat_i, c = highlight, maxr = 100)
    fig <- add_boxplot(fig, data = dfi, x = ~Var2, y = ~value, 
                       boxpoints = FALSE, color = ~Cat)
  }
  else if (nrow(mat_i) > 0) {
    if (nrow(mat_i) <= 3) {
      for (ii in seq_len(nrow(mat_i))) {
        fig <- add_lines(fig, x = colnames(mat_i), y = mat_i[ii, 
        ], name = rownames(mat_i)[ii], line = list(color = "rgb(150, 150, 150)", 
                                                   width = 1), showlegend = FALSE)
      }
    }
    dfi <- dfi <- convertDF(mat_i, c = highlight)
    fig <- add_markers(fig, data = dfi, x = ~Var2, y = ~value, 
                       color = ~Cat, text = paste0("<b>Feature: </b>", dfi$Var1, 
                                                   "<br>", "<b>Sample: </b>", tlpmap[dfi$Var2]), 
                       hoverinfo = "text")
  }
  fig <- plotly::layout(fig, xaxis = list(categoryarray = colnames(x), 
                                          categoryorder = "array", ticktext = o.name, tickvals = as.list(seq_len(ncol(x)) - 
                                                                                                           1), title = ""))
  if (is.null(extvar)) {
    fig <- plotly::layout(fig, yaxis = list(title = ylab))
    fig <- plotly::config(fig, toImageButtonOptions = list(format = "svg", 
                                                           filename = "omicsViewerPlot", width = 700, height = 700))
    return(fig)
  }
  figext <- plot_ly(x = colnames(x), y = extvar, type = "bar", 
                    showlegend = FALSE)
  ff <- subplot(figext, fig, shareX = TRUE, nrows = 2, heights = c(0.2, 
                                                                   0.8), margin = 0, titleY = TRUE)
  ff <- plotly::layout(ff, yaxis = list(title = ylab.extvar), 
                       yaxis2 = list(title = ylab))
  fig <- plotly::config(ff, toImageButtonOptions = list(format = "svg", 
                                                        filename = "omicsViewerPlot", width = 700, height = 700))
  toWebGL(fig)
}


plotly_boxplot_module <- function (input, output, session, reactive_param_plotly_boxplot, 
                                   reactive_checkpoint = reactive(TRUE)) 
{
  output$boxplotly <- renderPlotly({
    req(reactive_checkpoint())
    do.call(plotly_boxplot, args = reactive_param_plotly_boxplot())
  })
}


plotly_boxplot_ui <- function (id) 
{
  ns <- NS(id)
  shinycssloaders::withSpinner(plotlyOutput(ns("boxplotly")), 
                               type = 8, color = "green")
}


plotly_scatter <- function (x, y, xlab = "", ylab = "ylab", color = "", shape = "", 
                            size = 10, tooltips = NULL, regressionLine = FALSE, source = "scatterplotlysource", 
                            sizeRange = c(5, 15), highlight = NULL, highlightName = "Highlighted", 
                            inSelection = NA, vline = NULL, hline = NULL, rect = NULL, 
                            drawButtonId = NULL) 
{
  options(stringsAsFactors = FALSE)
  i1 <- (is.factor(x) || is.character(x) || is.logical(x)) && 
    is.numeric(y)
  i2 <- (is.factor(y) || is.character(y) || is.logical(y)) && 
    is.numeric(x)
  i3 <- is.numeric(y) && is.numeric(x)
  cutnumorchar <- function(x, n = 60, alt = "") {
    if (is.character(x) || is.factor(x)) {
      message("too many distinct values, not suitable for color mapping!")
      v <- alt
    }
    else if (is.numeric(x)) {
      v <- as.character(cut(x, breaks = n, include.lowest = TRUE, 
                            dig.lab = 3))
    }
    else stop("cutnumorchar: x needs to be one of objects: numeric, character, factor")
    v
  }
  if (length(unique(color)) > 60) 
    color <- cutnumorchar(color, alt = "")
  if (length(unique(shape)) > 60) 
    shape <- cutnumorchar(shape, alt = "")
  if (i1) {
    names(y) <- paste0("Y", seq_along(y))
    df <- beeswarm(y ~ x, corral = "wrap", do.plot = FALSE, 
                   corralWidth = 0.75)
    df$index <- fmatch(rownames(df), paste(x, names(y), sep = "."))
    tlp <- sprintf("<b>%s: </b>%s<br><b>%s: </b>%s", xlab, 
                   df$x.orig, ylab, signif(df$y, digits = 3))
    ixlab <- table(df$x.orig)
    ixlab <- structure(paste0(names(ixlab), " (n=", ixlab, 
                              ")"), names = names(ixlab))
    df$x.orig <- ixlab[df$x.orig]
  }
  else if (i2) {
    names(x) <- paste0("X", seq_along(x))
    df0 <- beeswarm(x ~ y, corral = "wrap", do.plot = FALSE, 
                    corralWidth = 0.9)
    df <- df0
    df$x <- df0$y
    df$y <- df0$x
    df$index <- fmatch(rownames(df), paste(y, names(x), sep = "."))
    tlp <- sprintf("<b>%s: </b>%s<br><b>%s: </b>%s", xlab, 
                   signif(df$x, digits = 3), ylab, df$y.orig)
  }
  else if (i3) {
    df <- data.frame(x = x, y = y)
    df$index <- seq_len(nrow(df))
    tlp <- sprintf("<b>%s: </b>%s<br><b>%s: </b>%s", xlab, 
                   signif(x, digits = 3), ylab, signif(y, digits = 3))
  }
  else {
    message("plotly_scatter: Unknown type x and/or y!")
    return(NULL)
  }
  f0 <- function(x, i) {
    if (length(x) == 1) 
      x <- rep(x, max(i, na.rm = TRUE))
    x[i]
  }
  df$color <- f0(color, i = df$index)
  df$shape <- f0(shape, i = df$index)
  df$size <- f0(size, i = df$index)
  if (!is.null(tooltips)) 
    tlp <- paste0("<b>", tooltips[df$index], "</b><br>", 
                  tlp)
  df$tlp <- tlp
  df$shape[is.na(df$shape)] <- "_NA_"
  df$color[is.na(df$color)] <- "_NA_"
  df$size[is.na(df$size)] <- min(df$size, na.rm = TRUE)
  if (!is.null(highlight)) {
    df$highlight <- FALSE
    df$highlight[which(df$index %in% highlight)] <- TRUE
  }
  df <- df[!(is.na(df$x) | is.na(df$y)), ]
  df <- df[order(df$x, decreasing = FALSE), ]
  df$xyid <- paste(df$x, df$y)
  if (nrow(df) == 0) 
    return(NULL)
  cc <- nColors(k = length(unique(df$color)))
  mop <- NA
  if (!is.na(inSelection)) {
    mop <- rep(0.2, nrow(df))
    mop[inSelection] <- 0.75
  }
  fig <- plot_ly(data = df, source = source)
  if (i1 || i2) {
    fig <- add_trace(fig, x = ~x, y = ~y, color = ~color, 
                     colors = cc, symbol = ~shape, size = ~size, sizes = sizeRange, 
                     marker = list(sizemode = "diameter", opacity = mop), 
                     type = "scatter", mode = "markers", text = ~tlp, 
                     hoverinfo = "text", showlegend = FALSE)
    if (i1) 
      fig <- plotly::layout(fig, yaxis = list(title = ylab), 
                            xaxis = list(title = xlab, tickvals = unique(round(df$x)), 
                                         ticktext = unique(df$x.orig)))
    if (i2) 
      fig <- plotly::layout(fig, xaxis = list(title = xlab), 
                            yaxis = list(title = ylab, tickvals = unique(round(df$y)), 
                                         ticktext = unique(df$x.orig)))
  }
  else {
    fig <- add_trace(fig, x = ~x, y = ~y, color = ~color, 
                     colors = cc, symbol = ~shape, size = ~size, sizes = sizeRange, 
                     marker = list(sizemode = "diameter", opacity = mop), 
                     text = ~tlp, hoverinfo = "text", type = "scatter", 
                     mode = "markers")
    if (regressionLine) {
      mod <- lm(df$y ~ df$x)
      prd <- predict(mod, newdata = data.frame(x), interval = "confidence")
      df <- cbind(df, prd)
      fig <- add_trace(fig, data = df, type = "scatter", 
                       mode = "lines", x = ~x, y = ~fit, line = list(color = "rgb(150, 150, 150)", 
                                                                     width = 1), hoverinfo = "skip", showlegend = FALSE, 
                       name = "regression<br>line")
      fig <- add_trace(fig, type = "scatter", mode = "lines", 
                       x = ~x, y = ~lwr, line = list(color = "transparent", 
                                                     width = 1), hoverinfo = "skip", showlegend = FALSE)
      fig <- add_trace(fig, type = "scatter", mode = "lines", 
                       x = ~x, y = ~upr, line = list(color = "transparent", 
                                                     width = 1), fill = "tonexty", fillcolor = "rgba(0,100,80,0.2)", 
                       hoverinfo = "skip", showlegend = FALSE)
      t <- list(size = 12)
      ct <- cor.test(x, y)
      p <- signif(ct$p.value, digits = 2)
      r <- signif(ct$estimate, digits = 2)
      fig <- plotly::layout(fig, title = list(text = sprintf("R = %s; p-value = %s", 
                                                             r, p), font = t))
    }
  }
  if (!is.null(highlight)) {
    fig <- add_trace(fig, x = df$x[df$highlight], y = df$y[df$highlight], 
                     type = "scatter", mode = "markers", name = highlightName, 
                     marker = list(size = 20, symbol = "circle-open", 
                                   color = "black", line = list(width = 3)), hoverinfo = "none", 
                     showlegend = FALSE)
  }
  if (is.list(rect)) {
    rect <- lapply(rect, function(r1) {
      if (any(!c("x0", "y0", "x1", "y1") %in% names(r1))) 
        return(NULL)
      list(type = "rect", fillcolor = "blue", line = list(color = "blue"), 
           opacity = 0.1, x0 = r1["x0"], x1 = r1["x1"], 
           xref = "x", y0 = r1["y0"], y1 = r1["y1"], yref = "y")
    })
  }
  fig <- plotly::layout(fig, xaxis = list(title = xlab), yaxis = list(title = ylab), 
                        shapes = rect, legend = list(orientation = "h", xanchor = "left", 
                                                     x = 0, y = 1, yanchor = "bottom"))
  if (is.null(drawButtonId)) 
    modeBarAdd <- NULL
  else modeBarAdd <- list(drawButton(drawButtonId))
  fig <- plotly::config(fig, toImageButtonOptions = list(format = "svg", 
                                                         filename = "omicsViewerPlot", width = 700, height = 700), 
                        modeBarButtonsToAdd = modeBarAdd)
  return(list(fig = toWebGL(fig), data = df))
}


plotly_scatter_module <- function (input,
                                   output, 
                                   session, 
                                   reactive_param_plotly_scatter, 
                                   reactive_regLine = reactive(FALSE),
                                   reactive_checkpoint = reactive(TRUE), 
                                   htest_var1 = reactive(NULL),
                                  htest_var2 = reactive(NULL)
                                  ) 
{
  options(warn = -1)
  ns <- session$ns
  hm <- reactive({
    x <- reactive_param_plotly_scatter()$x
    y <- reactive_param_plotly_scatter()$y
    i1 <- (is.factor(x) || is.character(x)) && is.numeric(y)
    i2 <- (is.factor(y) || is.character(y)) && is.numeric(x)
    i3 <- is.numeric(y) && is.numeric(x)
    list(scatter = i3, beeswarm = i1 || i2, beeswarm.vertical = i1)
  })
  choices <- reactive({
    req(reactive_checkpoint())
    req(hm()$beeswarm)
    x <- reactive_param_plotly_scatter()$x
    y <- reactive_param_plotly_scatter()$y
    if (hm()$beeswarm.vertical) {
      v <- x
      num <- y
    }
    else {
      v <- y
      num <- x
    }
    list(group = sort(unique(v)), x = num, f = v)
  })
  output$htest <- renderUI({
    fluidRow(column(3, uiOutput(ns("uiGroup1"))), column(3, 
                                                         uiOutput(ns("uiGroup2"))), column(6, DT::dataTableOutput(ns("testResult"))))
  })
  htv1 <- reactiveVal()
  htv2 <- reactiveVal()
  observe({
    if (is.null(htest_var1())) 
      htv1(choices()$group[1])
    else htv1(htest_var1())
    if (is.null(htest_var2())) 
      htv2(choices()$group[2])
    else htv2(htest_var2())
  })
  output$uiGroup1 <- renderUI({
    req(reactive_checkpoint())
    selectInput(inputId = ns("group1"), "group 1", choices = choices()$group, 
                selected = htv1(), selectize = TRUE, width = "100%")
  })
  output$uiGroup2 <- renderUI({
    req(reactive_checkpoint())
    selectInput(inputId = ns("group2"), "group 2", choices = choices()$group, 
                selected = htv2(), selectize = TRUE, width = "100%")
  })
  output$testResult <- DT::renderDataTable({
    req(reactive_checkpoint())
    req(input$group1)
    req(input$group2)
    req(x <- choices()$x)
    req(f <- choices()$f)
    r1 <- try(t.test(x[f == input$group1], x[f == input$group2]), 
              silent = TRUE)
    req(inherits(r1, "htest"))
    r2 <- wilcox.test(x[choices()$f == input$group1], x[choices()$f == 
                                                          input$group2])
    df <- data.frame(Diff = signif(r1$estimate[1] - r1$estimate[2], 
                                   digits = 3), `P t-test` = signif(r1$p.value, digits = 3), 
                     `P MVU-test` = signif(r2$p.value, digits = 3), check.names = FALSE, 
                     row.names = NULL)
    DT::datatable(df, options = list(searching = FALSE, lengthChange = FALSE, 
                                     dom = "t"), rownames = FALSE, class = "compact")
  })
  output$regTickBox <- renderUI({
    req(reactive_checkpoint())
    req(hm()$scatter)
    checkboxInput(ns("showRegLine"), "Regression line", reactive_regLine())
  })
  reactive_param_plotly_scatter_src <- reactive({
    x <- reactive_param_plotly_scatter()
    src <- ifelse(!is.null(x$source), x$source, "scatterplotly")
    x$source <- src
    x
  })
  plotter <- reactive({
    req(reactive_checkpoint())
    if (!hm()$scatter) {
      return(do.call(plotly_scatter, args = c(reactive_param_plotly_scatter_src(), 
                                              regressionLine = FALSE)))
    }
    req(!is.null(input$showRegLine))
    do.call(plotly_scatter, args = c(reactive_param_plotly_scatter_src(), 
                                     regressionLine = input$showRegLine))
  })
  output$plotly.scatter.output <- renderPlotly({
    req(plotter()$fig)
    plotter()$fig
  })
  rr <- reactive({
    id <- plotter()$data$index
    selected <- event_data("plotly_selected", source = reactive_param_plotly_scatter_src()$source)
    selected <- sort(id[fastmatch::"%fin%"(plotter()$data$xyid, 
                                           paste(selected$x, selected$y))])
    clicked <- event_data("plotly_click", source = reactive_param_plotly_scatter_src()$source)
    clicked <- sort(id[fastmatch::"%fin%"(plotter()$data$xyid, 
                                          paste(clicked$x, clicked$y))])
    list(selected = selected, clicked = clicked)
  })
  reactive({
    list(selected = rr()$selected, clicked = rr()$clicked, 
         regline = input$showRegLine, htest_V1 = input$group1, 
         htest_V2 = input$group2)
  })
}


plotly_scatter_ui <- function (id, height = "400px") 
{
  ns <- NS(id)
  tagList(uiOutput(ns("htest")), uiOutput(ns("regTickBox")), 
          shinycssloaders::withSpinner(plotlyOutput(ns("plotly.scatter.output"), 
                                                    height = height), type = 8, color = "green"))
}


prepOmicsViewer <- function (expr, pData, fData, PCA = TRUE, ncomp = min(8, ncol(expr)), 
                             pca.fillNA = TRUE, method = 'custom', t.test = NULL, ttest.fillNA = FALSE, ..., 
                             gs = NULL, stringDB = NULL, surv = NULL, SummarizedExperiment = TRUE) 
{
  print('Running Preomics')
  p0 <- pData
  de <- dim(expr)
  if (nrow(pData) != de[2]) 
    stop("nrow pData != ncol(exprs)")
  if (nrow(fData) != de[1]) 
    stop("nrow pData != nrow(exprs)")
  expr_rn <- rownames(expr)
  if (is.null(expr_rn) && !is.null(rownames(fData))) 
    expr_rn <- rownames(fData)
  if (is.null(expr_rn)) 
    expr_rn <- str_pad(seq_len(nrow(expr)), width = nchar(nrow(expr)), 
                       pad = "0")
  expr_rn <- make.names(expr_rn, unique = TRUE)
  expr_cn <- colnames(expr)
  if (is.null(expr_cn) && !is.null(rownames(pData))) 
    expr_cn <- rownames(pData)
  if (is.null(expr_cn)) 
    expr_cn <- str_pad(seq_len(ncol(expr)), width = nchar(ncol(expr)), 
                       pad = "0")
  expr_cn <- make.names(expr_cn, unique = TRUE)
  pData_rn <- rownames(pData)
  fData_rn <- rownames(fData)
  if (is.null(expr_rn) || is.null(expr_cn)) 
    stop("expr needs to be a matrix with dimnames")
  if (!identical(expr_cn, pData_rn)) {
    message("rownames of pData reset by colnames of exprs")
    rownames(pData) <- expr_cn
  }
  if (!identical(expr_rn, fData_rn)) {
    message("rownames of fData reset by rownames of exprs")
    rownames(fData) <- expr_rn
  }
  rownames(expr) <- expr_rn
  colnames(expr) <- expr_cn
  pData <- cbind(pData, numberOfFeatures = colSums(!is.na(expr)))
  colnames(pData) <- paste0("General|All|", trimws(colnames(pData)))
  colnames(fData) <- paste0("General|All|", trimws(colnames(fData)))
  if (PCA) {
    pc <- exprspca(expr, n = ncomp, fillNA = pca.fillNA, method = method)
    pData <- cbind(pData, pc$samples)
    fData <- cbind(fData, pc$features)
  }
  if (!is.null(t.test)) {
    print('Running for t-test')
    tres <- multi.t.test(x = expr, pheno = p0, compare = t.test, 
                         fillNA = ttest.fillNA, method = method, ...)
    fData <- cbind(fData, tres)
  }
  rk <- data.frame(apply(expr, 2, rank), stringsAsFactors = FALSE)
  colnames(rk) <- paste0("Rank|All|", colnames(rk))
  fData <- cbind(fData, rk)
  if (!is.null(stringDB)) {
    if (length(stringDB) == nrow(expr)) 
      strdb <- data.frame(`StringDB|All|ID` = stringDB, 
                          stringsAsFactors = FALSE, check.names = FALSE)
    else stop("stringDB should be a vector same length as nrow(expr), containing the IDs can be used to query stringDB")
    fData <- cbind(fData, strdb)
  }
  if (!is.null(surv)) {
    if (is.vector(surv) && length(surv) == nrow(expr)) {
      surv <- data.frame(`Surv|all|surv` = surv, stringsAsFactors = FALSE, 
                         check.names = FALSE)
    }
    else if (is.matrix(surv) || is.data.frame(surv)) {
      if (nrow(surv) != ncol(expr)) 
        stop("nrow(surv) should equals ncol(expr)")
      if (is.null(colnames(surv))) 
        colnames(surv) <- paste0("Surv", seq_len(ncol(surv)))
      surv <- data.frame(surv, stringsAsFactors = FALSE, 
                         check.names = TRUE)
      colnames(surv) <- paste0("Surv|all|", trimws(colnames(surv)))
    }
    else stop("incompatible 'surv'")
    sv <- vapply(surv, function(x) {
      all(grepl("\\D|+", sub("+$", "", surv)))
    }, logical(1))
    if (any(!sv)) 
      stop("only numbers and + are allowed in surv")
    pData <- cbind(pData, surv)
  }
  if (!is.null(gs)) {
    if (all(colnames(gs) %in% c("featureId", "gsId", "weight"))) {
      gs$featureId <- factor(rownames(fData)[gs$featureId])
      gs$gsId <- factor(gs$gsId)
      gs$weight <- as.integer(gs$weight)
    }
    else {
      if (is.vector(gs) && length(gs) == nrow(expr)) {
        gs <- matrix(gs, ncol = 1)
        colnames(gs) <- "geneset"
      }
      if (!inherits(gs, c("matrix", "dgCMatrix"))) 
        stop("gs should either be a (sparse) matrix or data.frame with three columns: featureId, gsId, weight!")
      if (is.null(rownames(gs))) 
        rownames(gs) <- rownames(fData)
      if (!(nrow(gs) == nrow(expr))) 
        stop("incompatible 'gs'")
      if (is.null(colnames(gs))) 
        stop("colnames of gs should not be null!")
      if (any(duplicated(colnames(gs)))) 
        stop("colnames of gs should be unique!")
    }
    attr(fData, "GS") <- gs
  }
  fx1 <- grep("ttest\\|(.*?)_vs_(.*?)\\|mean.diff", colnames(fData), 
              value = TRUE)
  fy1 <- intersect(colnames(fData), sub("mean.diff$", "log.fdr", 
                                        fx1))
  fx2 <- grep("PCA\\|All\\|PC1\\(", colnames(fData), value = TRUE)
  fy2 <- grep("PCA\\|All\\|PC2\\(", colnames(fData), value = TRUE)
  px <- grep("PCA\\|All\\|PC1\\(", colnames(pData), value = TRUE)
  py <- grep("PCA\\|All\\|PC2\\(", colnames(pData), value = TRUE)

  exprsWithAttr <- function(x, fillNA = FALSE, environment = FALSE, 
                            attrs = c("rowDendrogram", "colDendrogram")) {
    print(paste0(' #### using the imputation method from exprsWithAttr #### ',method))
    if (environment) 
      aenv <- new.env()
    else aenv <- list()
    mx <- as.matrix(x)
    for (i in attrs) attr(mx, i) <- attr(x, i)
    aenv$exprs <- mx
    if (fillNA) {
      mxf <- fillNA(mx, method = method)
      for (i in attrs) attr(mxf, i) <- attr(x, i)
      aenv$exprs_impute <- mxf
    }
    aenv
  }

  if (!SummarizedExperiment) {
    aenv <- exprsWithAttr(expr, fillNA = pca.fillNA || ttest.fillNA, 
                          environment = TRUE)
    res <- ExpressionSet(assayData = aenv, phenoData = AnnotatedDataFrame(pData), 
                         featureData = AnnotatedDataFrame(fData))
  }
  else {

    DataFrameWithAttr <- function(x) {
      attrs <- setdiff(names(attributes(x)), c("names", 
                                               "class", "row.names"))
      attr_list <- lapply(attrs, function(attr_name) attr(x, 
                                                          attr_name))
      names(attr_list) <- attrs
      x <- DataFrame(x, check.names = FALSE)
      for (i in names(attr_list)) attr(x, i) <- attr_list[[i]]
      x
    }


    aenv <- exprsWithAttr(expr, fillNA = pca.fillNA || ttest.fillNA, 
                          environment = FALSE)
    pd <- DataFrameWithAttr(pData)
    fd <- DataFrameWithAttr(fData)
    res <- SummarizedExperiment(assays = aenv, rowData = fd, 
                                colData = pd)
  }
  if (length(fx1) >= 1 && length(fy1) >= 1) {
    attr(res, "fx") <- fx1[1]
    attr(res, "fy") <- fy1[1]
  }
  else if (length(fx2) >= 1 && length(fy2) >= 1) {
    attr(res, "fx") <- fx2[1]
    attr(res, "fy") <- fy2[1]
  }
  if (length(px) >= 1 && length(py) >= 1) {
    attr(res, "sx") <- px[1]
    attr(res, "sy") <- py[1]
  }
  res
}


ptmotif_module <- function (input, output, session, pdata, fdata, expr, feature_selected, 
                            sample_selected, background) 
{
  ns <- session$ns
  triset <- reactive({
    req(fdata())
    i <- grep("^SeqLogo\\|", colnames(fdata()), value = TRUE)
    req(length(i) > 0)
    str_split_fixed(i, "\\|", n = 3)
  })
  xax <- reactiveVal()
  observe({
    xax(list(v1 = triset()[1, 1], v2 = triset()[1, 2], v3 = triset()[1, 
                                                                     3]))
  })
  v1 <- callModule(triselector_module, id = "tris_seqlogo", 
                   reactive_x = triset, label = "Sequence", reactive_selector1 = reactive(xax()$v1), 
                   reactive_selector2 = reactive(xax()$v2), reactive_selector3 = reactive(xax()$v3))
  scc <- reactive({
    req(v1())
    cs <- do.call(paste, list(v1(), collapse = "|"))
    req(cs %in% colnames(fdata()))
    cs
  })
  cleanSeqs <- function(x) {
    x <- unique(unlist(strsplit(x, ";")))
    x[which(nchar(x) > 0)]
  }
  bg.seqs <- reactiveVal(NULL)
  errText <- reactiveVal(NULL)
  observe({
    req(fdata())
    req(scc())
    bg.seqs(cleanSeqs(fdata()[, scc()]))
  })
  observe({
    req(bg.seqs())
    if (length(unique(nchar(bg.seqs()))) > 1) 
      errText("The length of sequences is different. The input of seqLogo analysis requires the sequences have the same length")
  })
  output$errorMsg <- renderText({
    if (is.null(errText())) 
      return(NULL)
    errText()
  })
  output$msg_ui <- renderUI({
    if (!is.null(errText())) 
      verbatimTextOutput(ns("errorMsg"))
    else tabsetPanel(tabPanel("Ratio selected/all", plotOutput(ns("plt"), 
                                                               height = "300px"), tabsetPanel(tabPanel("Selected seqs", 
                                                                                                       dataTableDownload_ui(ns("seqtable"))), tabPanel("Position weighted matrix", 
                                                                                                                                                       dataTableDownload_ui(ns("seqtable_rat"))))), tabPanel("Selected", 
                                                                                                                                                                                                             plotOutput(ns("plt.fg")), dataTableDownload_ui(ns("seqtable_fg"))), 
                     tabPanel("All", plotOutput(ns("plt.bg")), dataTableDownload_ui(ns("seqtable_bg"))))
  })
  foregroundSeqs <- reactive({
    req(fdata())
    req(scc())
    req(feature_selected())
    req(is.null(errText()))
    cleanSeqs(fdata()[feature_selected(), scc()])
  })
  bg.pfm <- reactive({
    req(bg.seqs())
    aaFreq(bg.seqs())
  })
  fg.pfm <- reactive({
    req(foregroundSeqs())
    aaFreq(foregroundSeqs())
  })
  logo <- reactive({
    req(bg.pfm())
    req(fg.pfm())
    motifRF(fg.pfm = fg.pfm(), bg.pfm = bg.pfm())
  })
  output$plt <- renderPlot({
    req(logo())
    ggseqlogo::ggseqlogo(data = logo()) + geom_vline(xintercept = (ncol(logo()) + 
                                                                     1)/2, linetype = "dashed", color = "orange", size = 1.5)
  })
  output$plt.fg <- renderPlot({
    req(d <- fg.pfm())
    ggseqlogo::ggseqlogo(data = d) + geom_vline(xintercept = (ncol(d) + 
                                                                1)/2, linetype = "dashed", color = "orange", size = 1.5)
  })
  output$plt.bg <- renderPlot({
    req(d <- bg.pfm())
    ggseqlogo::ggseqlogo(data = d) + geom_vline(xintercept = (ncol(d) + 
                                                                1)/2, linetype = "dashed", color = "orange", size = 1.5)
  })
  mat2df <- function(x) {
    data.frame(Name = rownames(x), x, stringsAsFactors = FALSE)
  }
  callModule(dataTableDownload_module, id = "seqtable", reactive_table = reactive({
    fg <- foregroundSeqs()
    fg <- fg[which(nchar(fg) > 0)]
    do.call(rbind, strsplit(fg, "|"))
  }), prefix = "motif", pageLength = 10)
  callModule(dataTableDownload_module, id = "seqtable_fg", 
             reactive_table = reactive(mat2df(fg.pfm())), prefix = "seqLogoPFM_foreground", 
             pageLength = 10)
  callModule(dataTableDownload_module, id = "seqtable_bg", 
             reactive_table = reactive(mat2df(bg.pfm())), prefix = "seqLogoPFM_background", 
             pageLength = 10)
  callModule(dataTableDownload_module, id = "seqtable_rat", 
             reactive_table = reactive(mat2df(logo())), prefix = "seqLogoPFM_ratio", 
             pageLength = 10)
}


ptmotif_ui <- function (id) 
{
  ns <- NS(id)
  tagList(triselector_ui(ns("tris_seqlogo")), uiOutput(ns("msg_ui")))
}


read_gmt <- function (x, id = NA, data.frame = FALSE) 
{
  x <- readLines(x)
  x <- strsplit(x, "\t")
  names <- vapply(x, function(xx) xx[1], character(1))
  x <- lapply(x, function(xx) {
    structure(xx[-(seq_len(2))], name = xx[1], link = xx[2])
  })
  names(x) <- names
  if (!is.na(id)) {
    id <- match.arg(id, c("SYMBOL", "ENTREZ"))
  }
  else {
    if (any(grepl("[A-Z]", unlist(x[seq_len(10)])))) 
      id <- "SYMBOL"
    else id <- "ENTREZ"
  }
  attr(x, "ID") <- id
  if (data.frame) {
    x <- data.frame(id = unlist(x), term = rep(names(x), 
                                               vapply(x, length, integer(1))), stringsAsFactors = FALSE)
  }
  x
}


read.proteinGroups <- function (x, quant = c("LF", "TMT")[1]) 
{
  func <- read.proteinGroups.lf
  getName <- function(x) {
    nm <- colnames(x$iBAQ)
    if (is.null(nm)) 
      nm <- colnames(x[[grep("LFQ", names(x))]])
    if (is.null(nm)) 
      nm <- colnames(x$Intensity)
    nm
  }
  if (quant == "TMT") {
    getName <- function(x) {
      colnames(v$Reporter.intensity.corrected)
    }
    func <- read.proteinGroups.tmt
  }
  v <- func(x)
  attr(v, "label") <- getName(v)
  id <- str_split_fixed(v$annot$Protein.IDs, pattern = ";", 
                        2)[, 1]
  gn <- str_split_fixed(v$annot$Gene.names, pattern = ";", 
                        2)[, 1]
  rn <- make.names(paste(gn, id, sep = "_"))
  ss <- names(which(unlist(vapply(v, nrow, FUN.VALUE = integer(1))) == 
                      nrow(v$annot)))
  for (i in ss) rownames(v[[i]]) <- rn
  v
}


read.proteinGroups.lf <- function (file) 
{
  pg <- read.delim(file, stringsAsFactors = FALSE)
  df <- data.frame(val = c("iBAQ.", "LFQ.intensity.", "Peptides.", 
                           "Razor...unique.peptides.", "Unique.peptides.", "Sequence.coverage.", 
                           "Intensity.", "MS.MS.Count.", "MS.MS.count.", "Identification.type."), 
                   log = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, 
                           FALSE, FALSE, FALSE), stringsAsFactors = FALSE)
  vi <- vapply(df$val, function(x) length(grep(x, colnames(pg))) > 
                 0, FUN.VALUE = logical(1))
  df <- df[vi, ]
  if (!"Intensity." %in% df$val) 
    stop("The proteinGroup.txt table should have at least Intensity columns!")
  i <- !(grepl("^REV_", pg$Majority.protein.IDs) | grepl("^CON_", 
                                                         pg$Majority.protein.IDs) | pg$Only.identified.by.site == 
           "+")
  annot <- pg[i, -grep(paste(df$val, collapse = "|"), colnames(pg))]
  getExpr <- function(x, type = "iBAQ.", log = TRUE, keep.row = NULL) {
    ic <- grep(type, colnames(x), ignore.case = FALSE, value = TRUE)
    ic <- setdiff(ic, "iBAQ.peptides")
    val <- apply(pg[, ic], 2, as.numeric)
    if (log) {
      val <- log10(val)
      val[is.infinite(val)] <- NA
    }
    if (!is.null(keep.row)) 
      val <- val[keep.row, ]
    colnames(val) <- gsub(type, "", colnames(val))
    val
  }
  ml <- mapply(function(val, log) getExpr(pg, type = val, log = log, 
                                          keep.row = i), val = df$val, log = df$log)
  names(ml) <- gsub("\\.$", "", df$val)
  ml$annot <- annot
  cnames <- colnames(ml$Intensity)
  for (i in names(ml)) {
    tmp <- ml[[i]]
    if (all(colnames(tmp) %in% cnames) && all(cnames %in% 
                                              colnames(tmp))) 
      ml[[i]] <- tmp[, cnames]
  }
  if (!is.null(ml$iBAQ)) {
    i <- which(rowSums(ml$iBAQ, na.rm = TRUE) == 0)
    if (length(i) > 0) 
      ml <- lapply(ml, function(x) x[-i, ])
    ml$iBAQ_mc <- sweep(ml$iBAQ, 2, matrixStats::colMedians(ml$iBAQ, 
                                                            na.rm = TRUE), "-") + median(ml$iBAQ, na.rm = TRUE)
  }
  ml
}


read.proteinGroups.tmt <- function (file, xref = NULL) 
{
  ab <- read.delim(file, stringsAsFactors = FALSE)
  ir <- c(grep("^CON_", ab$Majority.protein.IDs), grep("^REV_", 
                                                       ab$Majority.protein.IDs), which(ab$Only.identified.by.site == 
                                                                                         "+"))
  eSum <- c("Fraction", "Reporter.intensity.corrected", "Reporter.intensity", 
            "Reporter.intensity.count")
  ls <- list()
  for (i in eSum) {
    gb <- grep(paste0(i, ".[0-9]*$"), colnames(ab), value = TRUE)
    ls[[i]] <- apply(ab[-ir, gb, drop = FALSE], 2, as.numeric)
    ab[gb] <- NULL
  }
  lsind <- list()
  eInd <- c("Reporter.intensity.corrected", "Reporter.intensity.count", 
            "Reporter.intensity")
  for (i in eInd) {
    gb <- grep(i, colnames(ab), value = TRUE)
    lsind[[i]] <- apply(ab[-ir, gb, drop = FALSE], 2, as.numeric)
    ab[gb] <- NULL
  }
  lsind$Reporter.intensity.corrected.log10 <- log10(lsind$Reporter.intensity.corrected)
  lsind$Reporter.intensity.corrected.log10[is.infinite(lsind$Reporter.intensity.corrected.log10)] <- NA
  lsind$annot <- ab[-ir, ]
  lsind$Summed <- ls
  colnames(lsind$Reporter.intensity.corrected.log10) <- make.names(sub("Reporter.intensity.corrected.", 
                                                                       "", colnames(lsind$Reporter.intensity.corrected.log10)))
  ec <- intersect(eSum, names(lsind))
  for (i in ec) {
    nn <- sub(i, "", colnames(lsind[[i]]))
    nn <- make.names(sub("^.", "", nn))
    colnames(lsind[[i]]) <- nn
  }
  fn <- make.names(xref$label)
  if (!is.null(xref)) {
    lab <- make.names(paste(xref$channel, xref$mix))
    if (!identical(lab, colnames(lsind[[ec[[1]]]]))) 
      stop("columne does not match!")
    for (ii in ec) colnames(lsind[[ii]]) <- fn
    if (!is.null(lsind$Reporter.intensity.corrected.log10)) 
      colnames(lsind$Reporter.intensity.corrected.log10) <- fn
    lsind$xref <- xref
  }
  lsind
}


readESVObj  <- function (x) 
{
  if (grepl(".RDS$", x, ignore.case = TRUE)) {
    x <- asEsetWithAttr(readRDS(x))
    x <- tallGS(x)
  }
  else if (grepl("(.db|.sqlite|.sqlite3)$", x, ignore.case = TRUE)) {
    x <- dbConnect(RSQLite::SQLite(), x)
  }
  else stop("readESVObj: Unkown input format!")
  x
}


removeVarQC <- function (x, ref, positive = TRUE, ...) 
{
  ls <- list(...)
  if (length(ls) > 0) 
    x <- normalize.nQuantiles(x, ...)
  x0 <- x[, ref]
  decomp0 <- svd(x0)
  m0 <- decomp0$u %*% t(t(x) %*% decomp0$u)
  mm <- x - m0
  mm <- mm - rowMedians(mm)
  mm <- mm + rowMedians(x)
  if (positive) 
    mm[which(mm < 0)] <- 0
  mm
}


rowshift <- function (x, batch, ref = NULL, useMean = FALSE) 
{
  if (is.data.frame(x)) 
    x <- apply(x, 2, as.numeric)
  if (is.null(ref)) 
    ref <- seq_len(ncol(x))
  b_ref <- batch[ref]
  expr_ref <- x[, ref, drop = FALSE]
  grandmeans <- rowMeans(expr_ref, na.rm = TRUE)
  grandmeans2 <- rowMeans(x, na.rm = TRUE)
  for (i in unique(batch)) {
    off <- rowMeans(expr_ref[, b_ref == i, drop = FALSE], 
                    na.rm = TRUE) - grandmeans
    if (useMean) {
      off2 <- rowMeans(x[, batch == i, drop = FALSE], na.rm = TRUE) - 
        grandmeans2
      ona <- is.na(off)
      off[ona] <- off2[ona]
    }
    x[, batch == i] <- x[, batch == i] - off
  }
  x
}


sample_general_module <- function (input, output, session, reactive_phenoData, reactive_j = reactive(NULL), 
                                   reactive_status = reactive(NULL)) 
{
  ns <- session$ns
  triset <- reactive({
    trisetter(meta = reactive_phenoData(), combine = "none")
  })
  xax <- reactiveVal()
  v1 <- callModule(triselector_module, id = "tris_sample_general", 
                   reactive_x = triset, label = "Link selection to", reactive_selector1 = reactive(xax()$v1), 
                   reactive_selector2 = reactive(xax()$v2), reactive_selector3 = reactive(xax()$v3))
  attr4select_status <- reactiveVal()
  attr4select <- callModule(attr4selector_module, id = "a4_gp", 
                            reactive_meta = reactive_phenoData, reactive_triset = triset, 
                            reactive_status = attr4select_status)
  pheno <- reactive({
    req(v1()$variable)
    req(!v1()$variable %in% c("", "Select a variable!"))
    req(reactive_phenoData())
    cs <- do.call(paste, list(v1(), collapse = "|"))
    if (!cs %in% colnames(reactive_phenoData())) 
      return(NULL)
    val <- reactive_phenoData()[, cs]
    if (v1()$analysis == "Surv") {
      type <- "surv"
    }
    else if (is.numeric(val)) {
      type <- "beeswarm"
    }
    else if (is.character(val) || is.factor(val)) {
      type <- "table"
    }
    else {
      warnings("Unknown type of val: sample_general_module, return NULL!")
      return(NULL)
    }
    list(value = val, type = type)
  })
  select <- reactive({
    req(reactive_j())
    req(reactive_phenoData())
    select <- rep("Unselected", nrow(reactive_phenoData()))
    if (!is.null(reactive_j())) 
      select[rownames(reactive_phenoData()) %in% reactive_j()] <- "selected"
    select
  })
  output$sample_general_plot <- renderUI({
    req(pheno()$type)
    if (pheno()$type == "beeswarm") 
      r <- plotly_scatter_ui(ns("sample_general_beeswarm"))
    if (pheno()$type == "table") 
      r <- factorIndependency_ui(ns("sample_general_contab"))
    if (pheno()$type == "surv") 
      r <- survival_ui(ns("sample_general_surv"))
    tagList(column(11, r), column(1, attr4selector_ui(ns("a4_gp"), 
                                                      circle = FALSE, right = TRUE)))
  })
  htestV1 <- reactiveVal()
  htestV2 <- reactiveVal()
  vs_scatter <- callModule(plotly_scatter_module, id = "sample_general_beeswarm", 
                           reactive_param_plotly_scatter = reactive({
                             req(reactive_j())
                             req(pheno()$value)
                             tooltips <- attr4select$tooltips
                             if (is.null(tooltips)) 
                               tooltips <- rownames(reactive_phenoData())
                             l <- list(x = select(), y = pheno()$value, xlab = "", 
                                       ylab = do.call(paste, list(v1(), collapse = "|")), 
                                       tooltips = tooltips)
                             l$color <- attr4select$color
                             l$shape <- attr4select$shape
                             l$size <- attr4select$size
                             l$highlight <- attr4select$highlight
                             l$highlightName <- attr4select$highlightName
                             l
                           }), reactive_regLine = reactive(FALSE), reactive_checkpoint = reactive(pheno()$type == 
                                                                                                    "beeswarm"), htest_var1 = htestV1, htest_var2 = htestV2)
  callModule(factorIndependency_module, id = "sample_general_contab", 
             x = select, y = reactive(pheno()$value), reactive_checkpoint = reactive(pheno()$type == 
                                                                                       "table"))
  callModule(survival_module, id = "sample_general_surv", reactive_resp = reactive(pheno()$value), 
             reactive_strata = select, reactive_checkpoint = reactive(pheno()$type == 
                                                                        "surv"))
  metatab <- reactive({
    req(reactive_j())
    tab <- reactive_phenoData()
    tab <- tab[, grep("^General\\|", colnames(tab)), drop = FALSE]
    tab <- tab[reactive_j(), , drop = FALSE]
    ic <- vapply(tab, is.numeric, logical(1)) & vapply(tab, 
                                                       is.integer, logical(1))
    # tab[ic] <- lapply(tab[ic], signif, digits = 2)
    colnames(tab) <- sub("General\\|All\\|", "", colnames(tab))
    tab
  })
  callModule(dataTableDownload_module, id = "msatab", reactive_table = metatab, 
             prefix = "SampleTable_")
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status())) 
      return()
    xax(NULL)
    xax(list(v1 = s$xax[[1]], v2 = s$xax[[2]], v3 = s$xax[[3]]))
  })
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status())) 
      return()
    attr4select_status(NULL)
    attr4select_status(s$attr4)
  })
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status())) 
      return()
    htestV1(s$htestV1)
    htestV2(s$htestV2)
  })
  rv <- reactiveValues()
  observe(rv$xax <- v1())
  observe(rv$attr4 <- attr4select$status)
  observe({
    rv$htestV1 <- vs_scatter()$htestV1
    rv$htestV2 <- vs_scatter()$htestV2
  })
  reactive(reactiveValuesToList(rv))
}


sample_general_ui <- function (id) 
{
  ns <- NS(id)
  tagList(fluidRow(column(12, style = "margin-top: 0px;", triselector_ui(ns("tris_sample_general"), 
                                                                         right_margin = "5")), uiOutput(ns("sample_general_plot"))), 
          dataTableDownload_ui(ns("msatab")))
}

#nonstandardGenericFunction for "saveOmicsViewerDb" defined from package "omicsViewer"
saveOmicsViewerDb <-  function (obj, db.file, overwrite = TRUE) 
{
  standardGeneric("saveOmicsViewerDb")
}


shinyPlotTooltips <- function (input, output, session, points) 
{
  output$hover_info <- renderText({
    req(points())
    session$sendCustomMessage(type = "placeDraggable", message = list())
    paste0("<p style=\"background-color:rgb(250,250,250); padding:5px; border-radius:3px; border: 2px solid #CCC\">", 
           points(), "</p>")
  })
}


shinyPlotTooltipsUI  <-  function (id) 
{
  ns <- NS(id)
  htmlCode <- sprintf("// Get mouse coordinates\n    var mouseX, mouseY;\n    $(document).mousemove(function(e) {\n    mouseX = e.clientX;\n    mouseY = e.clientY;\n    }).mouseover();\n    \n    // Function to possition draggable, place on current mouse coordinates\n    Shiny.addCustomMessageHandler ('placeDraggable',function (message) {\n    var element = $('#%s').parent();\n    element.css({'top': mouseY + 'px', 'left' : mouseX + 'px'})\n    });", 
                      ns("hover_info"))
  tagList(tags$head(tags$script(HTML(htmlCode))), absolutePanel(fixed = TRUE, 
                                                                draggable = TRUE, htmlOutput(ns("hover_info"), style = "z-index: 9999;")))
}


sideCorKey  <-  function (x, label) 
{
  par(mar = c(0.1, 0.1, 2, 0.1))
  bp <- barplot(rep(1, length(x$key)), col = x$key, space = 0, 
                xaxs = "i", yaxs = "i", horiz = TRUE, border = x$key, 
                axes = FALSE, main = label)
  text(0.5, y = bp, names(x$key))
}


str2hclust  <-  function (x) 
{
  x <- strsplit(x, split = "__elementSplitter__")[[1]]
  ll <- list(merge = apply(read.delim(file = textConnection(x[1]), 
                                      stringsAsFactors = FALSE, header = FALSE), 2, as.integer), 
             height = as.numeric(strsplit(x[2], split = ";")[[1]]), 
             order = as.integer(strsplit(x[3], split = ";")[[1]]), 
             labels = strsplit(x[4], split = "=;=")[[1]], call = list("convertedHclustObj", 
                                                                      "str2hclust"), method = x[5], dist.method = x[6])
  if (length(ll$labels) == 0) 
    ll$labels <- NULL
  class(ll) <- "hclust"
  ll
}


string_module  <-  function (input, output, session, reactive_ids, reactive_status = reactive(NULL), 
                             active = reactive(FALSE)) 
{
  ns <- session$ns
  overflow <- reactive({
    length(reactive_ids()) > 300
  })
  output$error.msg <- renderText({
    sprintf("%s features selected [MAX 300 FEATURES ALLOWED!]", 
            length(reactive_ids()))
  })
  nk <- reactiveVal()
  observeEvent(input$run, {
    r <- stringNetwork(genes = reactive_ids(), taxid = input$tax)
    if (is.data.frame(r)) {
      if (nrow(r) > 999) {
        r <- r[order(r$score, decreasing = TRUE), ]
        r <- r[seq_len(999), ]
      }
    }
    nk(r)
  })
  gs <- reactiveVal()
  gs <- eventReactive(input$run, {
    req(!overflow())
    show_modal_spinner(text = "Querying database ...")
    tab <- stringGSA(genes = reactive_ids(), taxid = input$tax)
    remove_modal_spinner()
    if (inherits(tab, "character")) {
      return(tab)
    }
    colnames(tab) <- c("category", "term", "gene number", 
                       "background number", "TaxonId", "inputGenes", "preferredNames", 
                       "p value", "fdr", "description")
    tab
  })
  nores <- reactive({
    !is.data.frame(nk()) || !is.data.frame(gs())
  })
  output$nores.msg <- renderText({
    nores()
    req(nores())
    c(gs(), nk())[c(!is.data.frame(gs()), !is.data.frame(nk()))]
  })
  output$noresRet <- renderUI({
    verbatimTextOutput(ns("nores.msg"))
  })
  highlightP <- reactiveVal(1)
  output$network <- renderForceNetwork({
    req(!nores())
    req(nrow(nk()) > 0)
    stringD3Net(ntwk = nk(), gsa = gs(), i = highlightP(), 
                label = input$showLabel)
  })
  tt <- callModule(dataTableDownload_module, id = "strtab", 
                   reactive_table = eventReactive(gs(), {
                     req(!nores())
                     req(!overflow())
                     req(nrow(gs()) > 0)
                     gs()[, c("category", "term", "gene number", "background number", 
                              "p value", "fdr", "description")]
                   }), prefix = "FeatureTable_")
  observe({
    req(tt())
    highlightP(tt())
  })
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status())) 
      return()
    updateTextInputIcon(session, "tax", value = s$tax)
    updateCheckboxInput(session, "showLabel", value = s$showLabel)
    if (active()) 
      shinyjs::click("run")
  })
  reactive({
    list(tax = input$tax, showLabel = input$showLabel)
  })
}


string_ui  <-  function (id) 
{
  ns <- NS(id)
  tagList(fluidRow(column(4, offset = 0, style = "padding-left:15px; padding-right:2px; padding-top:0px; padding-bottom:0px", 
                          textInputIcon(inputId = ns("tax"), label = NULL, value = "9606", 
                                        icon = list("Taxomony Code"))), column(6, offset = 0, 
                                                                               style = "padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px", 
                                                                               verbatimTextOutput(ns("error.msg"))), column(2, offset = 0, 
                                                                                                                            style = "padding-left:15px; padding-right:2px; padding-top:0px; padding-bottom:0px", 
                                                                                                                            actionButton(ns("run"), "Run!"))), uiOutput(ns("noresRet")), 
          dataTableDownload_ui(ns("strtab")), checkboxInput(ns("showLabel"), 
                                                            label = "Show labels", value = FALSE), forceNetworkOutput(ns("network")))
}


stringD3Net  <-  function (ntwk, gsa, i, label = FALSE) 
{
  nd <- data.frame(name = unique(unlist(ntwk[, c("preferredName_A", 
                                                 "preferredName_B")])), stringsAsFactors = FALSE)
  rownames(nd) <- nd$name
  nd$group <- 1
  i <- strsplit(gsa$preferredNames[i], ",")[[1]]
  nd$group[nd$name %in% i] <- 2
  links <- data.frame(source = fmatch(ntwk$preferredName_A, 
                                      nd$name) - 1, target = fmatch(ntwk$preferredName_B, nd$name) - 
                        1, value = (ntwk$score - 0.4)^2 * 10)
  colorfunc <- networkD3::JS("colorfunc = function(i) { return i == 1 ? \"#64A0C8\" : \"#E37222\" };")
  lab <- ifelse(label, 1, 0)
  forceNetwork(Links = links, Nodes = nd, Source = "source", 
               linkColour = "gray", Target = "target", Value = "value", 
               NodeID = "name", charge = -5, Group = "group", colourScale = colorfunc, 
               opacity = 0.7, fontSize = 12, opacityNoHover = lab, zoom = TRUE, 
               legend = TRUE)
}


stringGSA  <-  function (genes, taxid = 9606, background = NULL, backgroundStringId = FALSE, 
                         caller = "omicsViewer") 
{
  string_api_url <- "https://string-db.org/api"
  output_format <- "tsv"
  method <- "enrichment"
  request_url <- paste(string_api_url, output_format, method, 
                       sep = "/")
  params <- list(species = taxid, caller_identity = caller)
  if (!is.null(background)) {
    if (!backgroundStringId) {
      g <- getStringId(genes, taxid = taxid)
      if (inherits(g, "character")) 
        return(g)
      sg <- g$stringId
      genes <- intersect(genes, g$queryItem)
    }
    else {
      sg <- background
      genes <- intersect(genes, background)
    }
    if (length(genes) == 0) 
      stop("No genes exist in the current background!")
    params$background_string_identifiers <- paste(g$stringId, 
                                                  collapse = "%0d")
  }
  params$identifiers <- paste(genes, collapse = "%0d")
  response <- httr::GET(url = request_url, query = params)
  dd <- httr::content(response)
  if (!is.data.frame(dd)) {
    message("stringdb does not return valid results, return NULL!")
    return(dd)
  }
  as.data.frame(dd)
}


stringNetwork  <-  function (genes, taxid = 9606, caller = "omicsViewer") 
{
  string_api_url <- "https://string-db.org/api"
  output_format <- "tsv"
  method <- "network"
  if (inherits(genes, "list")) 
    genes <- na.omit(unlist(genes))
  genes <- gsub("#|%", "", genes)
  request_url <- paste(string_api_url, output_format, method, 
                       sep = "/")
  params <- list(identifiers = paste(genes, collapse = "%0d"), 
                 species = taxid, required_score = 500, caller_identity = caller)
  params <- paste(vapply(names(params), function(x) paste(x, 
                                                          params[[x]], sep = "="), FUN.VALUE = character(1)), collapse = "&")
  results <- GET(paste(request_url, params, sep = "?"))
  dd <- httr::content(results)
  dd <- unique(dd)
  if (!is.data.frame(dd)) {
    message("stringdb does not return valid results, return NULL!")
    return(dd)
  }
  as.data.frame(dd)
}


survival_module  <-  function (input, output, session, reactive_resp, reactive_strata, 
                               reactive_checkpoint = reactive(TRUE)) 
{
  ns <- session$ns
  dat <- reactive({
    req(reactive_checkpoint())
    y <- reactive_resp()
    data.frame(time = as.numeric(sub("\\+$", "", y)), event = as.integer(grepl("\\+$", 
                                                                               y)), strata = reactive_strata(), stringsAsFactors = FALSE)
  })
  output$censor_output <- renderUI({
    nm <- max(dat()$time, na.rm = TRUE)
    fluidRow(column(12, offset = 0, style = "padding-left:5px; padding-right:5px; padding-top:0px; padding-bottom:0px", 
                    div(style = "display: inline-block;vertical-align:top;", 
                        h5("Censor at")), div(style = "padding-left:25px; display: inline-block;vertical-align:top; width:65%;", 
                                              sliderInput(ns("censor"), label = NULL, min = min(dat()$time, 
                                                                                                na.rm = TRUE), max = nm, value = nm))))
  })
  output$kmplot <- renderPlot({
    req(input$censor)
    df <- dat()
    i <- which(df$time > input$censor)
    df$time[i] <- input$censor
    df$event[i] <- 0
    fit <- survfit(Surv(time, event) ~ strata, data = df)
    lab <- ""
    if (length(df$strata > 1)) {
      r <- surv_pvalue(fit, data = df, method = "survdiff")
      lab <- paste(r$method, r$pval.txt)
    }
    ggsurvplot(fit, data = df, risk.table = TRUE, conf.int = TRUE, 
               pval = lab, surv.median.line = "hv")
  })
}


survival_ui <- function (id) 
{
  ns <- NS(id)
  tagList(uiOutput(ns("censor_output")), plotOutput(ns("kmplot")))
}

tallGS <- function (obj) 
{
  fd <- Biobase::fData(obj)
  scn <- str_split_fixed(colnames(fd), "\\|", n = 3)
  ir <- which(scn[, 1] == "GS")
  gscsc <- attr(fd, "GS")
  if (length(ir) > 0) {
    igs <- fd[, ir, drop = FALSE]
    colnames(igs) <- make.names(scn[ir, 3])
    gs <- totall(igs)
    fd <- fd[, -ir]
    attr(fd, "GS") <- gs
  }
  else if (inherits(gscsc, c("dgCMatrix", "lgCMatrix"))) {
    gs <- csc2list(gscsc)
    attr(fd, "GS") <- gs
  }
  Biobase::fData(obj) <- fd
  obj
}


text2num <- function (x) 
{
  if (is.null(x)) 
    return(NULL)
  x0 <- try(eval(parse(text = x)), silent = TRUE)
  if (!is.numeric(x0)) {
    warning("text2num: cannot convert x to num!")
    return(NULL)
  }
  x0
}



totall <- function (gsmat) 
{
  gs <- as.matrix(gsmat)
  gs[gs == 0] <- NA
  gs <- reshape2::melt(gs, na.rm = TRUE)
  colnames(gs) <- c("featureId", "gsId", "weight")
  gs
}

triselector_module <- function (input, output, session, reactive_x, reactive_selector1 = reactive(NULL), 
                                reactive_selector2 = reactive(NULL), reactive_selector3 = reactive(NULL), 
                                label = "Group Label:") 
{
  ns <- session$ns
  output$groupLabel <- renderUI({
    h5(HTML(sprintf("<b>%s</b>", label)))
  })
  inte <- reactive({
    aa <- grepl("tris_feature_general", ns("x"))
    length(aa) > 0 && aa
  })
  observeEvent(list(reactive_selector1()), {
    req(reactive_x())
    if (length(names(input)) == 0) 
      return(NULL)
    cc <- unique(reactive_x()[, 1])
    if (!is.null(reactive_selector1())) 
      ss <- reactive_selector1()
    else ss <- cc[1]
    updateSelectInput(session, inputId = "analysis", choices = cc, 
                      selected = ss)
  })
  observeEvent(list(names(input), reactive_x()), {
    req(reactive_x())
    if (length(names(input)) == 0) 
      return(NULL)
    cc <- unique(reactive_x()[, 1])
    if (input$analysis %in% cc) 
      ss <- input$analysis
    else if (!is.null(reactive_selector1())) 
      ss <- reactive_selector1()
    else ss <- cc[1]
    updateSelectInput(session, inputId = "analysis", choices = cc, 
                      selected = ss)
  })
  observeEvent(input$analysis, {
    req(reactive_x())
    req(input$analysis == "")
    cc <- unique(reactive_x()[, 1])
    if (is.null(reactive_selector1())) 
      ss <- cc[1]
    else ss <- reactive_selector1()
    updateSelectInput(session, inputId = "analysis", choices = cc, 
                      selected = ss)
  })
  observe({
    req(reactive_x())
    input$analysis
    req(input$analysis)
    cc <- unique(reactive_x()[reactive_x()[, 1] == input$analysis, 
                              2])
    updateSelectInput(session, inputId = "subset", choices = cc, 
                      selected = reactive_selector2())
  })
  observe({
    input$analysis
    input$subset
    req(input$analysis)
    req(input$subset)
    req(reactive_x())
    cc <- reactive_x()[, 3][reactive_x()[, 1] == input$analysis & 
                              reactive_x()[, 2] == input$subset]
    cc <- c("Select a variable!", cc)
    preselected <- try(match.arg(reactive_selector3(), cc), 
                       silent = TRUE)
    if (inherits(preselected, "try-error")) 
      preselected <- NULL
    updateSelectInput(session, inputId = "variable", choices = cc, 
                      selected = preselected)
  })
  reactive({
    req(input$variable)
    list(analysis = input$analysis, subset = input$subset, 
         variable = input$variable)
  })
}


triselector_ui <- function (id, right_margin = "20") 
{
  rmar <- sprintf("padding-left:2px; padding-right:%spx; padding-top:2px; padding-bottom:2px", 
                  right_margin)
  ns <- NS(id)
  tagList(fluidRow(column(2, offset = 0, align = "right", style = "padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px", 
                          uiOutput(ns("groupLabel"))), column(3, offset = 0, style = "padding:2px;", 
                                                              selectInput(inputId = ns("analysis"), label = NULL, choices = NULL, 
                                                                          selectize = TRUE, width = "100%")), column(4, offset = 0, 
                                                                                                                     style = "padding:2px;", selectInput(inputId = ns("subset"), 
                                                                                                                                                         label = NULL, choices = NULL, selectize = TRUE, width = "100%")), 
                   column(3, offset = 0, style = rmar, selectInput(inputId = ns("variable"), 
                                                                   label = NULL, choices = NULL, selectize = TRUE, width = "100%"))))
}


trisetter <- function (meta, expr = NULL, combine) 
{
  req(meta)
  nm <- colnames(meta)
  cgs <- attr(meta, "GS")
  if (!is.null(cgs)) 
    cgs <- paste("GS", "All", unique(cgs$gsId), sep = "|")
  if (!is.null(expr)) {
    combine <- combine[1]
    combine <- match.arg(combine, choices = c("pheno", "feature", 
                                              "none"))
    if (combine == "feature") {
      if (nrow(expr) != nrow(meta)) 
        stop("nrow(meta) needs to equal nrow(expr) when combined by 'feature'!")
      rt <- paste0("Sample|Auto|", colnames(expr))
    }
    else if (combine == "pheno") {
      if (ncol(expr) != nrow(meta)) 
        stop("nrow(meta) needs to equal ncol(expr) when combined by 'pheno'!")
      rt <- paste0("Feature|Auto|", rownames(expr))
    }
    else {
      rt <- NULL
    }
    nm <- c(nm, cgs, rt)
  }
  v <- str_split_fixed(nm, "\\|", n = 3)
  attr(v, "seed") <- Sys.time()
  v
}


validMQFolder <- function (dir) 
{
  l <- list(valid = FALSE)
  i1 <- file.exists(file.path(dir, "mqpar.xml"))
  i2 <- list.files(dir, pattern = "^proteinGroups.txt$", recursive = TRUE, 
                   full.names = TRUE)
  if (i1 & length(i2) == 1) {
    l$valid <- TRUE
    l$mqpar <- file.path(dir, "mqpar.xml")
    l$txt <- dirname(i2)
    l$filePath <- file.path(l$txt, "proteinGroups.txt")
  }
  l
}


value2color <- function (x, n = 10) 
{
  if (is.numeric(x)) {
    nl <- as.factor(x)
    if (length(unique(x)) > n) 
      nl <- cut(x, n)
    cp <- colorRampPalette(brewer.pal(n = 7, name = "Blues"))(n)
    names(cp) <- levels(nl)
    cc <- cp[nl]
  }
  else if (is.factor(x) || is.character(x) || is.logical(x)) {
    x <- as.character(x)
    nl <- sort(unique(x))
    if (length(nl) > 40) {
      cp <- nColors(k = 40)
      cp <- rep(cp, ceiling(length(nl)/40))[seq_along(nl)]
    }
    else {
      cp <- nColors(k = length(nl))
    }
    names(cp) <- nl
    cp <- cp[sort(names(cp))]
    cc <- cp[x]
  }
  else stop("x is not one of: 1) numeric vectoer 2) char vector 3) factor 4) logical!")
  cc[is.na(cc)] <- "white"
  list(color = cc, key = cp)
}


varSelector <- function (x, expr, meta, alternative = NULL) 
{
  if (is.null(x)) 
    return(NULL)
  if (x$variable %in% c("Select a variable!", "")) {
    if (is.null(alternative)) 
      return(NULL)
    if (alternative %in% c("", "Select a variable!")) 
      return(NULL)
    return(alternative)
  }
  lab <- paste(x, collapse = "|")
  if (x$analysis == "GS") {
    gs <- attr(meta, "GS")
    if (is.null(gs)) 
      return(NULL)
    id1 <- gs[which(gs$gsId == x$variable), ]
    if (nrow(id1) == 0) 
      return(NULL)
    vec <- rep(0, nrow(meta))
    vec[fmatch(id1$featureId, rownames(meta))] <- id1$weight
    return(vec)
  }
  selectIfCan <- function(x, i, dim) {
    if (dim == 1) {
      if (!all(i %in% rownames(x))) 
        return(NULL)
      r <- x[i, ]
    }
    else if (dim == 2) {
      if (!all(i %in% colnames(x))) 
        return(NULL)
      r <- x[, i]
    }
    else stop("dim should be either 1 or 2!")
    r
  }
  if (x$analysis == "Feature" && x$subset == "Auto") {
    x <- selectIfCan(expr, x$variable, dim = 1)
  }
  else if (x$analysis == "Sample" && x$subset == "Auto") {
    x <- selectIfCan(expr, x$variable, dim = 2)
  }
  else {
    x <- selectIfCan(meta, lab, dim = 2)
  }
  if (!is.null(x)) 
    attr(x, "label") <- lab
  x
}


vectORA <- function (gs, i, background = NA, minOverlap = 3, minSize = 5, 
                     maxSize = Inf, gs_desc = NULL, feature_desc = NULL, unconditional.or = TRUE, 
                     mtc.method = "fdr", sort = c("none", "p.value", "OR")[1]) 
{
  if (is.na(background)) 
    background <- nrow(gs)
  mat <- gs
  mat[gs != 0] <- 1
  mat[is.na(mat)] <- 0
  ncand <- colSums(mat != 0)
  ic <- which(ncand >= minSize & ncand <= maxSize)
  if (length(ic) == 0) {
    message("No gene set in range c(minSize, maxSize), return NULL!")
    return(NULL)
  }
  mat <- mat[, ic, drop = FALSE]
  nol <- colSums(mat[i, ])
  ic <- which(nol >= minOverlap)
  if (length(ic) == 0) {
    message("No gene set have overlap >= minOverlap, return NULL!")
    return(NULL)
  }
  mat <- mat[, ic, drop = FALSE]
  nol <- nol[ic]
  bdf <- vectORA.core(n.overlap = nol, n.de = rep(length(i), 
                                                  ncol(mat)), n.gs = colSums(mat), n.bkg = background, 
                      unconditional.or = unconditional.or, mtc.method = mtc.method)
  gs_annot_x <- ""
  if (!is.null(gs_desc)) 
    gs_annot_x <- gs_desc[colnames(mat)]
  if (is.null(feature_desc)) 
    feature_desc <- rownames(gs)
  if (is.matrix(mat)) {
    fi <- feature_desc[i]
    overlap <- lapply(seq_len(ncol(mat)), function(j) {
      x <- mat[i, j, drop = FALSE]
      fi[x != 0]
    })
  }
  else if (inherits(mat, "dgCMatrix")) {
    overlap <- csc2list(mat[i, ])
    overlap <- split(as.character(overlap$featureId), overlap$gsId)
  }
  else stop("vectORA should be either a matrix or dgCMatrix!")
  rs <- cbind(pathway = colnames(mat), desc = gs_annot_x, bdf)
  rs$overlap_ids <- overlap
  sort <- sort[1]
  if (sort == "p.value") {
    rs <- rs[order(rs$p.value, decreasing = FALSE), ]
  }
  else if (sort == "OR") {
    rs <- rs[order(rs$OR, decreasing = TRUE), ]
  }
  else if (sort != "none") 
    warning("Unknown sort method, the results are not sorted!")
  rs
}


vectORA.core <-function (n.overlap, n.de, n.gs, n.bkg, unconditional.or = TRUE, 
                         mtc.method = "fdr") 
{
  pval <- phyper(q = n.overlap - 1, m = n.gs, n = n.bkg - n.gs, 
                 k = n.de, lower.tail = FALSE)
  if (length(pval) == 0) 
    return(data.frame(p.value = numeric(0), p.adjusted = numeric(0), 
                      OR = numeric(0), size_overlap = numeric(0), size_geneset = numeric(0), 
                      size_input = numeric(0), size_backgroung = numeric(0), 
                      stringsAsFactors = FALSE))
  if (unconditional.or) 
    or <- (n.overlap/(n.de - n.overlap))/((n.gs - n.overlap)/(n.bkg - 
                                                                n.gs - n.de + n.overlap))
  else {
    or <- function(n.overlap, n.gs, n.de, n.bkg) {
      m <- n.gs
      n <- n.bkg - n.gs
      k <- n.de
      x <- n.overlap
      lo <- pmax(0L, k - n)
      hi <- pmin(k, m)
      supportl <- mapply(":", lo, hi, SIMPLIFY = FALSE)
      vapply(seq_len(length(x)), function(i) {
        support <- supportl[[i]]
        logdc <- dhyper(support, m[i], n[i], k[i], log = TRUE)
        dnhyper <- function(ncp) {
          d <- logdc + log(ncp) * support
          d <- exp(d - max(d))
          d/sum(d)
        }
        mnhyper <- function(ncp) {
          if (ncp == 0) 
            return(lo[i])
          if (ncp == Inf) 
            return(hi[i])
          sum(support * dnhyper(ncp))
        }
        mle <- function(x) {
          if (x == lo[i]) 
            return(0)
          if (x == hi[i]) 
            return(Inf)
          mu <- mnhyper(1)
          if (mu > x) 
            uniroot(function(t) mnhyper(t) - x, c(0, 
                                                  1))$root
          else if (mu < x) 
            1/uniroot(function(t) mnhyper(1/t) - x, c(.Machine$double.eps, 
                                                      1))$root
          else 1
        }
        mle(x[i])
      }, numeric(1))
    }
    or <- or(n.overlap = n.overlap, n.gs = n.gs, n.de = n.de, 
             n.bkg = n.bkg)
  }
  data.frame(p.value = pval, p.adjusted = p.adjust(pval, method = mtc.method), 
             OR = or, size_overlap = n.overlap, size_geneset = n.gs, 
             size_input = n.de, size_backgroung = n.bkg, stringsAsFactors = FALSE)
}


vectORATall <- function (gs, i, background, minOverlap = 2, minSize = 2, maxSize = Inf, 
                         gs_desc = NULL, feature_desc = NULL, unconditional.or = TRUE, 
                         mtc.method = "fdr", sort = c("none", "p.value", "OR")[1]) 
{
  if (!is.factor(gs$featureId) || !is.factor(gs$gsId)) 
    stop("gs should have two columns of FACTOR named as featureId and gsId!")
  cnt_gs <- table(gs$gsId)
  cnt_gs <- cnt_gs[cnt_gs >= minSize & cnt_gs <= maxSize]
  if (length(cnt_gs) == 0) {
    message("No gene set in range c(minSize, maxSize), return NULL!")
    return(NULL)
  }
  gs <- gs[gs$gsId %in% names(cnt_gs), ]
  gsi <- gs[gs$featureId %in% i, ]
  cnt_gsi <- table(gsi$gsId)
  cnt_gsi <- cnt_gsi[cnt_gsi >= minOverlap]
  if (length(cnt_gsi) == 0) {
    message("No gene set have overlap >= minOverlap, return NULL!")
    return(NULL)
  }
  gsi <- gsi[gsi$gsId %in% names(cnt_gsi), ]
  bdf <- vectORA.core(n.overlap = c(cnt_gsi), n.de = rep(length(i), 
                                                         length(cnt_gsi)), n.gs = c(cnt_gs[names(cnt_gsi)]), n.bkg = background, 
                      unconditional.or = unconditional.or, mtc.method = mtc.method)
  gs_annot_x <- ""
  if (!is.null(gs_desc)) 
    gs_annot_x <- gs_desc[names(cnt_gsi)]
  if (is.null(feature_desc)) 
    feature_desc <- as.character(gsi$featureId)
  rs <- cbind(pathway = names(cnt_gsi), desc = gs_annot_x, 
              bdf)
  overlap <- split(feature_desc, gsi$gsId)[as.character(rs$pathway)]
  rs$overlap_ids <- unname(overlap)
  sort <- sort[1]
  if (sort == "p.value") {
    rs <- rs[order(rs$p.value, decreasing = FALSE), ]
  }
  else if (sort == "OR") {
    rs <- rs[order(rs$OR, decreasing = TRUE), ]
  }
  else if (sort != "none") 
    warning("Unknown sort method, the results are not sorted!")
  rs
}



