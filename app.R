#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)

library(shinybusy)
library(tidyverse)
library(ggplot2)


library(corrplot)



# File size limit increased
options(shiny.maxRequestSize = 100*1024^2)



# Define UI for application that draws a histogram
ui <- dashboardPage(
  dashboardHeader(title = "GSE107509"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Calculator", tabName = "calculator", icon = icon("calculator")),
      menuItem("Principal Component Analysis", tabName = "pca", icon = icon("chart-area", lib = "font-awesome")),
      menuItem("Random Forest", tabName = "rf", icon = icon("project-diagram"))
    )
  ),
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "calculator",
              fluidRow(
                box(
                  width = 12,
                  status = "primary",
                  
                  h2("Gene Expression Data"),
                  
                  br(),
                  
                  fluidRow(
                    valueBox("GSE107509", "MICROARRAY DATA", icon = icon("database"), color = "light-blue"),
                    valueBox(659, "SAMPLES", icon = icon("user"), color = "light-blue"),
                    valueBox(54715, "GENES", icon = icon("dna"), color = "light-blue"),
                    valueBox("GENE EXPRESSION", "DESCRIPTION", icon = icon("th"), color = "light-blue"),
                    valueBox(166, "OUTCOME: SUBCLINICAL ACUTE REJECTION", icon = icon("times"), color = "light-blue"),
                    valueBox(493, "OUTCOME: TRANSPLANT EXCELLENCE", icon = icon("check"), color = "light-blue")
                  ),
                )
              ),
              
              
              fluidRow(
                box(
                  width = 12,
                  status = "primary",
                  
                  h2("Top 12 Genes (By Variance)"),
                  
                  br(),
                  
                  # Sample selector
                  h4("Select a Sample:"),
                  sliderInput("sample_selection", label = NULL,
                              min = 1, max = 659, value = 1),
                  
                  valueBoxOutput("top1"),  # --------------
                  valueBoxOutput("top2"),  # --------------
                  valueBoxOutput("top3"),  # --------------
                  valueBoxOutput("top4"),  # --------------
                  valueBoxOutput("top5"),  # --------------
                  valueBoxOutput("top6"),  # --------------
                  valueBoxOutput("top7"),  # --------------
                  valueBoxOutput("top8"),  # --------------
                  valueBoxOutput("top9"),  # --------------
                  valueBoxOutput("top10"),  # --------------
                  valueBoxOutput("top11"),  # --------------
                  valueBoxOutput("top12"),  # --------------
                )
              ),
              
              fluidRow(
                box(
                  width = 12,
                  status = "primary",
                  
                  h2("Outcome"),
                  
                  br(),
                  
                  valueBoxOutput("sample_outcome"),  # --------------
                )
              )
      ), # End Calculator Tab
      
      tabItem(tabName = "pca",
              fluidRow(
                box(
                  width = 12,
                  status = "primary",
                  
                  h2("Principal Component Analysis"),
                  
                  br(),
                  
                  
                  plotOutput("pca_variance_explained_plot"),
                  sliderInput("select_pca_number", "Number of PCA:",
                              min = 1, max = 659,
                              value = 50),
                  valueBoxOutput("pca_num_box"),
                  valueBoxOutput("pca_value_box")
                )
              )
      ), # End PCA tab
      
      tabItem(tabName = "rf",  
              
              fluidRow(
                box(
                  width = 12,
                  status = "primary",
                  
                  h2("Random Forest Evaluation"),
                  
                  br(),
                  
                  # Value Boxes
                  valueBox("Random Forest", "Model", icon = icon("", lib = "font-awesome"), color = "light-blue"), 
                  valueBox(500, "Number of Trees", icon = icon("", lib = "font-awesome"), color = "light-blue"), 
                  valueBox("Accuracy", "Metric", icon = icon("", lib = "font-awesome"), color = "light-blue"),  
                  valueBox("Cross Validation", "Validation Technique", icon = icon("check", lib = "font-awesome"), color = "light-blue"),  
                  valueBox(5, "Number of Folds", icon = icon("map", lib = "font-awesome"), color = "light-blue"),  
                  valueBox(25, "Number of Repetitions", icon = icon("redo-alt", lib = "font-awesome"),color = "light-blue"), 
                  
                )
              ),
              fluidRow(
                box(
                  width = 12,
                  status = "primary",
                  
                  h2("How Does PCA Affect Random Forest Accuracy?"),
                  
                  br(),
                  
                  selectizeInput("pca_components_input", "Select Number of Principal Components :", choices = c("5", "10", "25", "50", "75", "100", "150", "200"),  multiple = TRUE,
                                 selected = c("5", "10"), options = list(placeholder = 'Select number of PCs to compare')),
                  plotOutput("pca_rf_accuracy_plot")
                )
              )
      ) # End RF tab
      
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$selected_pcas <- renderPrint({
    input$pca_components_input
    print(input$pca_components_input)
  })
  
  # Gene Expression Data for 659 samples
  GSE107509 = reactive({
    
    # Start progress bar
    progress <- shiny::Progress$new()
    progress$set(message = "Loading Gene Expression Data", value = 0)
    
    # Directory for GSE107509
    directory = "data/GSE107509_Split/"
    
    # Read in the files
    fileNames <- list.files(directory)
    
    # Read in all files to make a table
    GSE107509 = as.data.frame(readr::read_csv("data/GSE107509_Split/GSE107509_GeneNames.txt"))
    
    # Skip First (GeneName) data
    for(i in 2:length(fileNames)){
      print(fileNames[i])
      progress$inc(1/7, detail = "Concatenating")
      temptable <- readr::read_csv(file.path(directory, fileNames[i]))
      # Concatenate the second column (This particular person)
      GSE107509 <- cbind(GSE107509, temptable)
    }
    
    # Get gene names
    gene_names = GSE107509$Gene
    
    # Remove first column - The X column we saved
    GSE107509 = GSE107509[,-1]
    
    # Set gene names as rows
    rownames(GSE107509) = gene_names
    
    progress$inc(1, detail = "Done")
    
    return(GSE107509)
    
  })
  
  # Transposes Gene Expression Data
  GSE107509_t = reactive({
    progress <- shiny::Progress$new()
    progress$set(message = "Transposing GSE107509", value = 0)
    
    gse_transpose = as.data.frame(t(GSE107509()))
    
    progress$inc(1, detail = "Done")
    return(gse_transpose)
  })
  
  # Mean for each gene
  gse_mean = reactive({
    return(read.csv("data/GSE107509_mean.csv"))
  })
  
  # Var for each gene
  gse_var = reactive({
    return(read.csv("data/GSE107509_var.csv"))
  })
  
  # Rejection Status for each sample
  gse_rejection_status = reactive({
    data = read.csv("data/rejection_status_GSE107509.txt")
    data = data %>% select(outcome)
    return(data[input$sample_selection,])
  })
  
  # Top 12 gene variances
  top_12_var = reactive({
    return(gse_var() %>% arrange(desc(Var))  %>% head(12))
  })
  
  # Returns gene data for selected sample
  selected_sample <- reactive({
    data = data.frame(Gene = rownames(GSE107509()), Value = GSE107509()[,input$sample_selection])
    rownames(data) = rownames(GSE107509())
    return(data)
  })
  
  # Get PCA variance table
  pca_variance = reactive({
    return(read.csv("data/GSE107509_pca_variance_explained.csv"))
  })
  
  # PCA explained value from selection
  pca_explained_value = reactive({
    value = pca_variance() %>% filter(PCA_Number == input$select_pca_number)
    return(value$Variance_Explained[1])
  })
  
  
  # Random Forest Output
  precomputed_random_forest_accuracies = reactive({
    # Read in all files to make a table
    data = as.data.frame(readr::read_csv("data/GSE107509_precomputed_rf_acc.txt"))
    data = data %>% filter(PCA %in% input$pca_components_input)
    
    data$PCA = as.factor(data$PCA)
    
    print(str(data))
    return(data)
  })
  
  output$pca_rf_accuracy_plot = renderPlot({
    plot = ggplot(data = precomputed_random_forest_accuracies(), aes(x = PCA, y = Accuracy, fill = PCA)) + geom_boxplot() + theme_minimal() + labs (title = "Number of Principal Components vs Random Forest Accuracy", x = "PCA", y = "Accuracy") 

    return(plot)
  })
  
  output$pca_variance_explained_plot <- renderPlot({
    plot = ggplot(data = pca_variance(), aes(x = PCA_Number, y = Variance_Explained)) + geom_line() + theme_minimal() + geom_vline(xintercept=input$select_pca_number, color='coral') + 
      labs(title = "Cumulative Sum of Variance Explained for PCA", x = "Number of Principal Components", y = "Variance Explained (%)")
    return(plot)
  })
  
  
  
  output$pca_value_box <- renderValueBox({
    valueBox(
      paste0(round(pca_explained_value() * 100, 2), "%"), "Variance Explained", icon = icon("percent", lib = "font-awesome"),
      color = "green"
    )
  })
  
  output$pca_num_box <- renderValueBox({
    valueBox(
      input$select_pca_number, "Number of Principal Components", icon = icon("", lib = "font-awesome"),
      color = "green"
    )
  })
  
  
  
  # --- Top 12 Genes:
  
  output$top1 <- renderValueBox({
    
    # Get top gene variance (1)
    gene = top_12_var()$Gene[1]
    
    # Get the value of the gene expression for this sample
    this_gene_value = selected_sample() %>% filter(Gene == gene)
    this_gene_value = round(this_gene_value$Value[1],4)
    
    # Get the mean value gene expression
    mean_gene_value = round(top_12_var()$Var[1], 4)
    
    # Value box visualisation based on if this > mean or <
    if (this_gene_value > mean_gene_value) {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-up", lib = "font-awesome"),
        color = "green"
      )
    }
    else {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-down", lib = "font-awesome"),
        color = "red"
      )
    }
  })
  
  output$top2 <- renderValueBox({
    
    # Get top gene variance (2)
    gene = top_12_var()$Gene[2]
    
    # Get the value of the gene expression for this sample
    this_gene_value = selected_sample() %>% filter(Gene == gene)
    this_gene_value = round(this_gene_value$Value[1],4)
    
    # Get the mean value gene expression
    mean_gene_value = round(top_12_var()$Var[2], 4)
    
    # Value box visualisation based on if this > mean or <
    if (this_gene_value > mean_gene_value) {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-up", lib = "font-awesome"),
        color = "green"
      )
    }
    else {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-down", lib = "font-awesome"),
        color = "red"
      )
    }
  })
  
  output$top3 <- renderValueBox({
    
    # Get top gene variance (3)
    gene = top_12_var()$Gene[3]
    
    # Get the value of the gene expression for this sample
    this_gene_value = selected_sample() %>% filter(Gene == gene)
    this_gene_value = round(this_gene_value$Value[1],4)
    
    # Get the mean value gene expression
    mean_gene_value = round(top_12_var()$Var[3], 4)
    
    # Value box visualisation based on if this > mean or <
    if (this_gene_value > mean_gene_value) {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-up", lib = "font-awesome"),
        color = "green"
      )
    }
    else {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-down", lib = "font-awesome"),
        color = "red"
      )
    }
  })
  
  output$top4 <- renderValueBox({
    
    # Get top gene variance (4)
    gene = top_12_var()$Gene[4]
    
    # Get the value of the gene expression for this sample
    this_gene_value = selected_sample() %>% filter(Gene == gene)
    this_gene_value = round(this_gene_value$Value[1],4)
    
    # Get the mean value gene expression
    mean_gene_value = round(top_12_var()$Var[4], 4)
    
    # Value box visualisation based on if this > mean or <
    if (this_gene_value > mean_gene_value) {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-up", lib = "font-awesome"),
        color = "green"
      )
    }
    else {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-down", lib = "font-awesome"),
        color = "red"
      )
    }
  })
  
  output$top5 <- renderValueBox({
    
    # Get top gene variance (5)
    gene = top_12_var()$Gene[5]
    
    # Get the value of the gene expression for this sample
    this_gene_value = selected_sample() %>% filter(Gene == gene)
    this_gene_value = round(this_gene_value$Value[1],4)
    
    # Get the mean value gene expression
    mean_gene_value = round(top_12_var()$Var[5], 4)
    
    # Value box visualisation based on if this > mean or <
    if (this_gene_value > mean_gene_value) {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-up", lib = "font-awesome"),
        color = "green"
      )
    }
    else {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-down", lib = "font-awesome"),
        color = "red"
      )
    }
  })
  
  output$top6 <- renderValueBox({
    
    # Get top gene variance (6)
    gene = top_12_var()$Gene[6]
    
    # Get the value of the gene expression for this sample
    this_gene_value = selected_sample() %>% filter(Gene == gene)
    this_gene_value = round(this_gene_value$Value[1],4)
    
    # Get the mean value gene expression
    mean_gene_value = round(top_12_var()$Var[6], 4)
    
    # Value box visualisation based on if this > mean or <
    if (this_gene_value > mean_gene_value) {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-up", lib = "font-awesome"),
        color = "green"
      )
    }
    else {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-down", lib = "font-awesome"),
        color = "red"
      )
    }
  })
  
  output$top7 <- renderValueBox({
    
    # Get top gene variance (7)
    gene = top_12_var()$Gene[7]
    
    # Get the value of the gene expression for this sample
    this_gene_value = selected_sample() %>% filter(Gene == gene)
    this_gene_value = round(this_gene_value$Value[1],4)
    
    # Get the mean value gene expression
    mean_gene_value = round(top_12_var()$Var[7], 4)
    
    # Value box visualisation based on if this > mean or <
    if (this_gene_value > mean_gene_value) {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-up", lib = "font-awesome"),
        color = "green"
      )
    }
    else {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-down", lib = "font-awesome"),
        color = "red"
      )
    }
  })
  
  output$top8 <- renderValueBox({
    
    # Get top gene variance (8)
    gene = top_12_var()$Gene[8]
    
    # Get the value of the gene expression for this sample
    this_gene_value = selected_sample() %>% filter(Gene == gene)
    this_gene_value = round(this_gene_value$Value[1],4)
    
    # Get the mean value gene expression
    mean_gene_value = round(top_12_var()$Var[8], 4)
    
    # Value box visualisation based on if this > mean or <
    if (this_gene_value > mean_gene_value) {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-up", lib = "font-awesome"),
        color = "green"
      )
    }
    else {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-down", lib = "font-awesome"),
        color = "red"
      )
    }
  })
  
  output$top9 <- renderValueBox({
    
    # Get top gene variance (9)
    gene = top_12_var()$Gene[9]
    
    # Get the value of the gene expression for this sample
    this_gene_value = selected_sample() %>% filter(Gene == gene)
    this_gene_value = round(this_gene_value$Value[1],4)
    
    # Get the mean value gene expression
    mean_gene_value = round(top_12_var()$Var[9], 4)
    
    # Value box visualisation based on if this > mean or <
    if (this_gene_value > mean_gene_value) {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-up", lib = "font-awesome"),
        color = "green"
      )
    }
    else {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-down", lib = "font-awesome"),
        color = "red"
      )
    }
  })
  
  output$top10 <- renderValueBox({
    
    # Get top gene variance (10)
    gene = top_12_var()$Gene[10]
    
    # Get the value of the gene expression for this sample
    this_gene_value = selected_sample() %>% filter(Gene == gene)
    this_gene_value = round(this_gene_value$Value[1],4)
    
    # Get the mean value gene expression
    mean_gene_value = round(top_12_var()$Var[10], 4)
    
    # Value box visualisation based on if this > mean or <
    if (this_gene_value > mean_gene_value) {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-up", lib = "font-awesome"),
        color = "green"
      )
    }
    else {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-down", lib = "font-awesome"),
        color = "red"
      )
    }
  })
  
  output$top11 <- renderValueBox({
    
    # Get top gene variance (11)
    gene = top_12_var()$Gene[11]
    
    # Get the value of the gene expression for this sample
    this_gene_value = selected_sample() %>% filter(Gene == gene)
    this_gene_value = round(this_gene_value$Value[1],4)
    
    # Get the mean value gene expression
    mean_gene_value = round(top_12_var()$Var[11], 4)
    
    # Value box visualisation based on if this > mean or <
    if (this_gene_value > mean_gene_value) {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-up", lib = "font-awesome"),
        color = "green"
      )
    }
    else {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-down", lib = "font-awesome"),
        color = "red"
      )
    }
  })
  
  output$top12 <- renderValueBox({
    
    # Get top gene variance (12)
    gene = top_12_var()$Gene[12]
    
    # Get the value of the gene expression for this sample
    this_gene_value = selected_sample() %>% filter(Gene == gene)
    this_gene_value = round(this_gene_value$Value[1],4)
    
    # Get the mean value gene expression
    mean_gene_value = round(top_12_var()$Var[12], 4)
    
    # Value box visualisation based on if this > mean or <
    if (this_gene_value > mean_gene_value) {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-up", lib = "font-awesome"),
        color = "green"
      )
    }
    else {
      valueBox(
        this_gene_value, paste0("GENE: ", gene), icon = icon("chevron-down", lib = "font-awesome"),
        color = "red"
      )
    }
  })
  
  output$sample_outcome <- renderValueBox({
    
    # Get outcome from selected sample
    outcome = gse_rejection_status()
    
    # Value box visualisation based on outcome
    if (outcome == "TRANSPLANT EXCELLENCE") {
      valueBox(
        outcome, "OUTCOME", icon = icon("check", lib = "font-awesome"),
        color = "green"
      )
    }
    else {
      valueBox(
        "REJECTION", "OUTCOME", icon = icon("times", lib = "font-awesome"),
        color = "red"
      )
    }
  })
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
