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

library(corrplot)



# File size limit increased
options(shiny.maxRequestSize = 100*1024^2)



# Define UI for application that draws a histogram
ui <- dashboardPage(
    dashboardHeader(title = "Predictor"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Calculator", tabName = "calculator", icon = icon("calculator")),
            menuItem("Data Explorer", tabName = "data_explorer", icon = icon("project-diagram")),
            menuItem("Model Details", tabName = "model_details", icon = icon("project-diagram"))
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
                        )
                    )
                ),
                fluidRow(
                    box(
                        title = "INPUT",
                        width = 4,
                        status = "primary",
                        
                        
                        # Load Gene Data
                        actionButton(inputId = "load_gene_data", label = "Load Gene Data"),
                                
                        # Input: Select a file ----
                        fileInput("file1", "Upload Patient CSV File",
                                  multiple = TRUE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        
                        # Horizontal line ----
                        tags$hr(),
                        
                        # Input: Checkbox if file has header ----
                        checkboxInput("header", "Header", TRUE),
                        
                        # Input: Select separator ----
                        radioButtons("sep", "Separator",
                                     choices = c(Comma = ",",
                                                 Tab = "\t"),
                                     selected = "\t"),
                        
                        # Horizontal line ----
                        tags$hr(),
                        
                        h3("Data Preview"),
                        DT::dataTableOutput("contents") 
                                
                            
                    ),
                    box(
                        title = "OUTPUT",
                        width = 8,
                        status = "primary",
                        
                        h3("Top 6 Genes (By Variance)"),
                        DT::dataTableOutput("top_6_var_table"),
                        
                        valueBoxOutput("top1"),  # --------------
                        valueBoxOutput("top2"),  # --------------
                        valueBoxOutput("top3"),  # --------------
                        valueBoxOutput("top4"),  # --------------
                        valueBoxOutput("top5"),  # --------------
                        valueBoxOutput("top6"),  # --------------
                        
                        DT::dataTableOutput("GSE107509_Table"),   # --- Need to add this back in if want to use gse
                        DT::dataTableOutput("gse_table_rejection_status")
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
            tabItem(tabName = "model_details",
                    DT::dataTableOutput("gse_table_mean"),
                    DT::dataTableOutput("gse_table_var"),
                    DT::dataTableOutput("gse_table_mean_var")
                    
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
      return(read.csv("data/rejection_status_GSE107509.txt"))
    })
    
    # Top 6 gene variances
    top_6_var = reactive({
      return(gse_var() %>% arrange(desc(Var))  %>% head(6))
    })
    
    # Gets data from uploaded file
    uploaded_file = reactive({
      
      # input$file1 will be NULL initially. After the user selects
      # and uploads a file, head of that data file by default,
      # or all rows if selected, will be shown.
      req(input$file1)
      
      progress <- shiny::Progress$new()
      progress$set(message = "Loading Data", value = 0)
      
      df <- read.csv(input$file1$datapath,
                     header = input$header,
                     sep = input$sep)
      progress$inc(1, detail = "Done")
      #progress$close()
      
      return(df)
    })
    
    
    
    output$GSE107509_Table <- DT::renderDataTable({
      return(GSE107509())
    })
    
    output$top_6_var_table <- DT::renderDataTable({
      return(top_6_var())
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
    
    
    
    observeEvent(input$load_gene_data, {
        print("Clicked Button")
        gse()
    })

    output$contents <- DT::renderDataTable({
        print("Loading Selected File")
        return(uploaded_file())
    })
    
    output$gse_table <- DT::renderDataTable({
        print("Loading GSE File")
        return(gse())
    })
    
    
    # --- Top 6 Genes:
    
    get_top_1_gene = reactive({
      
      val = uploaded_file() %>% filter(Gene == "224588_PM_at")
      
      print(val)
      val = val$Value[1]
      
      print(val)
      
      return(val)
      
    })
    
    output$top1 <- renderValueBox({
      
      # Get top gene variance (1)
      gene = top_6_var()$Gene[1]
      
      # Get the value of the gene expression for this sample
      this_gene_value = uploaded_file() %>% filter(Gene == "224588_PM_at")
      this_gene_value = round(this_gene_value$Value[1],4)
      
      # Get the mean value gene expression
      mean_gene_value = round(top_6_var()$Var[1], 4)
      
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
      gene = top_6_var()$Gene[2]
      
      # Get the value of the gene expression for this sample
      this_gene_value = uploaded_file() %>% filter(Gene == gene)
      this_gene_value = round(this_gene_value$Value[1],4)
      
      # Get the mean value gene expression
      mean_gene_value = round(top_6_var()$Var[2], 4)
      
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
      gene = top_6_var()$Gene[3]
      
      # Get the value of the gene expression for this sample
      this_gene_value = uploaded_file() %>% filter(Gene == gene)
      this_gene_value = round(this_gene_value$Value[1],4)
      
      # Get the mean value gene expression
      mean_gene_value = round(top_6_var()$Var[3], 4)
      
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
      gene = top_6_var()$Gene[4]
      
      # Get the value of the gene expression for this sample
      this_gene_value = uploaded_file() %>% filter(Gene == gene)
      this_gene_value = round(this_gene_value$Value[1],4)
      
      # Get the mean value gene expression
      mean_gene_value = round(top_6_var()$Var[4], 4)
      
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
      gene = top_6_var()$Gene[5]
      
      # Get the value of the gene expression for this sample
      this_gene_value = uploaded_file() %>% filter(Gene == gene)
      this_gene_value = round(this_gene_value$Value[1],4)
      
      # Get the mean value gene expression
      mean_gene_value = round(top_6_var()$Var[5], 4)
      
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
      gene = top_6_var()$Gene[6]
      
      # Get the value of the gene expression for this sample
      this_gene_value = uploaded_file() %>% filter(Gene == gene)
      this_gene_value = round(this_gene_value$Value[1],4)
      
      # Get the mean value gene expression
      mean_gene_value = round(top_6_var()$Var[6], 4)
      
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
