library(shiny)
source("FoldChange_calculation.R")
# Define server logic for random distribution application
shinyServer(function(input, output) {

  output$mRNAplot <- renderPlot({    
    data = read.csv("mRNA_pathway-31Jan2014.csv", header = T)
    data1 = data[-1, -1]
    data = as.numeric(as.matrix(data1))
    data = matrix(data, nrow = nrow(data1), ncol = nrow(data1))
    image(t(data),
         main=paste('SIMS table'), xlab = "genes", ylab = "Patient ID")
  })
  
  # Generate a summary of the data
  datatable_mrna<-reactive({
    inputfile_mrna<-input$mytable
    if(is.null(inputfile_mrna)){return(NULL)}
    mrna_dtable<-read.csv(inputfile_mrna$datapath, sep=",", header= TRUE)    
    mrna_dtable = mrna_dtable #[-1,-1]
    return(mrna_dtable)
  })
  
  datatable_mirna <- reactive({
    inputfile_mirna <- input$mytable
    if(is.null(inputfile_mirna)){return(NULL)}
    mirna_dtable <- read.csv(inputfile_mirna$datapath, sep= ",", header = TRUE)
    return(mirna_dtable)
  })
  output$mRNAcontents<-renderTable({
    if (is.null(input$mytable)){return(NULL)}
    return(datatable_mrna()[, 1:5])
  })
  
  output$miRNAcontents <- renderTable({
    if(is.null(input$mytable)){return(NULL)}
    return(datatable_mirna()[, 1:5])
  })
  
  output$mRNAsummary<-renderPrint({
    if (is.null(input$mytable)){return("File upload error: No file uploaded")}
    final_data_gene_expression = calculate_fold_change(datatable_mrna(), "mRNA")
    return(final_data_gene_expression[1:5, 1:5])
  })
  
  output$miRNAsummary <- renderPrint({
    if(is.null(input$mytable)){return("File upload error: No file uploaded")}
    return(summary(datatable_mirna()[, 1:5]))
  })
    
    
    
    
    
    
    
})
