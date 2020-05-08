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
      menuItem("Data Explorer", tabName = "data_explorer", icon = icon("project-diagram")),
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
                  title = "GENE EXPRESSION DATA",
                  width = 12,
                  status = "primary",
                  
                  fluidRow(
                    valueBox("GSE107509", "MICROARRAY DATA", icon = icon("database"), color = "light-blue"),
                    valueBox(659, "SAMPLES", icon = icon("user"), color = "light-blue"),
                    valueBox(54715, "GENES", icon = icon("dna"), color = "light-blue"),
                    valueBox("GENE EXPRESSION", "DESCRIPTION", icon = icon(""), color = "light-blue"),
                    valueBox(166, "OUTCOME: SUBCLINICAL ACUTE REJECTION", icon = icon("times"), color = "light-blue"),
                    valueBox(493, "OUTCOME: TRANSPLANT EXCELLENCE", icon = icon("check"), color = "light-blue")
                  ),
                )
              ),
              
              
              fluidRow(
                box(
                  title = "TOP 12 GENES (BY VARIANCE)",
                  width = 12,
                  status = "primary",
                  
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
                  title = "OUTCOME",
                  width = 12,
                  status = "primary",
                  
                  valueBoxOutput("sample_outcome"),  # --------------
                )
              )
      ), # End Calculator Tab
      tabItem(tabName = "data_explorer",
              
              fluidRow(
                box(
                  title = "CORRELATIONS",
                  width = 12,
                  status = "primary",
                  
                  # Select Number of Correlations 
                  sliderInput("corr_number", "Number correlations (First n/54715):",
                              min = 0, max = 50, value = 5),
                  
                  plotOutput(outputId = "corrplot")
                ),
                
                box(
                  title = "BOXPLOT GENE EXPRESSION",
                  width = 12,
                  status = "primary",
                  # Issue with step? It is always 2x whatever it is set to 
                  numericInput("sample_number1", "Sample Number 1:", min = 1, max = 659, value = 1,  width = '25%'),
                  numericInput("sample_number2", "Sample Number 2:", min = 1, max = 659, value = 2, width = '25%'),
                  plotOutput(outputId = "gene_expression_comparison_boxplot")
                )
              )
              
      ),
      tabItem(tabName = "pca",
              fluidRow(
                box(
                  title = "PRINCIPAL COMPONENT ANALYSIS",
                  width = 12,
                  status = "primary",
                  
                  plotOutput("pca_variance_explained_plot"),
                  sliderInput("select_pca_number", "Number of PCA:",
                              min = 1, max = 659,
                              value = 50),
                  valueBoxOutput("pca_num_box"),
                  valueBoxOutput("pca_value_box")
                )
              )
      ),
      
      tabItem(tabName = "rf",  
      )
      
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
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
  
  
   
  
  
  # Compute Principal Component Analysis
  pca_reactive = reactive({
    
    progress <- shiny::Progress$new()
    progress$set(message = "Calculating PCA", value = 0)
    
    pca = prcomp(GSE107509_t())
    
    progress$inc(1, detail = "Done")
    
    return(pca)
  })
  
  output$pca <- renderText({
    print("Loading Selected File")
    return(pca_reactive())
  })
  
  pca_variance = reactive({
    return(read.csv("data/GSE107509_pca_variance_explained.csv"))
  })
  
  output$pca_variance_explained_plot <- renderPlot({
    plot = ggplot(data = pca_variance(), aes(x = PCA_Number, y = Variance_Explained)) + geom_line() + theme_minimal() + geom_vline(xintercept=input$select_pca_number, color='coral') + 
      labs(title = "Cumulative Sum of Variance Explained for PCA", x = "Number of Principal Components", y = "Variance Explained (%)")
    return(plot)
  })
  
  pca_explained_value = reactive({
    value = pca_variance() %>% filter(PCA_Number == input$select_pca_number)
    return(value$Variance_Explained[1])
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
  
  
  
  output$GSE107509_Table <- DT::renderDataTable({
    
    return(GSE107509())
  })
  
  output$top_6_var_table <- DT::renderDataTable({
    return(top_12_var())
  })
  
  output$gse_table_rejection_status <- DT::renderDataTable({
    return(gse_rejection_status())
  })
  
  
  output$gse_table_mean <- DT::renderDataTable({
    return(data.frame(gse_mean()))
  })
  
  output$gse_table_var <- DT::renderDataTable({
    return(data.frame(gse_var()))
  })
  
  output$gse_table_mean_var <- DT::renderDataTable({
    return(data.frame(gse_mean(), Var = gse_var()$Var))
  })

  
  
  
  
  output$contents <- DT::renderDataTable({
    print("Loading Selected File")
    return(uploaded_file())
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
        color = "green", width = 12
      )
    }
    else {
      valueBox(
        "REJECTION", "OUTCOME", icon = icon("times", lib = "font-awesome"),
        color = "red", width = 12
      )
    }
  })
  
  
  # --- 
  
  
  output$corrplot = renderPlot({
    cor_matrix<-cor(GSE107509_t()[1:input$corr_number])
    diag(cor_matrix)<-0
    return(corrplot(cor_matrix, method="square"))
  })
  
  
  output$gene_expression_comparison_boxplot = renderPlot({
    
    p <- ggplot(melt(GSE107509()[,c(input$sample_number1, input$sample_number2)]), aes(x=variable, y=value)) +  
      geom_boxplot(outlier.colour="black", outlier.shape=16,
                   outlier.size=0.5, notch=FALSE) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs (x = "Sample", y = "Expression Value") + theme_minimal() + coord_flip()
    return(p)
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
