library(shiny)
library(plyr)
#setwd("/srv/shiny-server/sims/")
#library(googleVis)
#require(graphics)
# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 30MB.
options(shiny.maxRequestSize = 30*1024^2)
source("FoldChange_calculation.R")
source("CopyNumberAlteration.R")
source("mutation_score.R")
# Define server logic for random distribution application
shinyServer(function(input, output) { 
  
  # Generate a summary of the data
  datatable_mrna<-reactive({
    inputfile_mrna<-input$mytable_mrna
    if(is.null(inputfile_mrna)){return(NULL)}
    mrna_dtable<-read.csv(inputfile_mrna$datapath, header = input$header, sep = input$sep, stringsAsFactors = F)    # header for header file, sep for separator
    return(mrna_dtable)
  })
  
  histology <- reactive({
    Histology <- c("AC","AC","SCC","AC","AC","SCC","AC","AC","SCC","SCC","LCC","AC","SCC","AC","SCC","AC","AC","AC","Other SCLC","LCC","AC","LCC","LCC","AC","AC","SCC","SCC","SCC","AC","SCC","AC","AC","Other ADEC","SCC","SCC","AC","SCC","LCC","AC","AC","SCC","AC","AC","AC","SCC","AC","SCC","SCC","SCC","AC","SCC","AC","SCC","AC","AC","SCC","AC","AC","AC","AC","SCC","SCC","AC","AC","SCC","SCC","AC","LCC","SCC","AC","AC","AC","LCC","SCC","SCC","SCC","SCC","AC","LCC","SCC","AC","SCC","AC","SCC","SCC","SCC","AC","SCC","LCC","SCC","AC","AC","AC","AC","LCC","AC","Other ADEC","AC","SCC","SCC","AC","LCC","SCC","AC","SCC","AC","SCC","LCC","AC","SCC","AC","SCC","AC","SCC","SCC","SCC","SCC","SCC","SCC","AC","AC")
    return(Histology)
  })
  
  datatable_mirna<-reactive({
    inputfile_mirna<-input$mytable_mirna
    if(is.null(inputfile_mirna)){return(NULL)}
    mirna_dtable<-read.csv(inputfile_mirna$datapath, header = input$header, sep = input$sep)    
    return(mirna_dtable)
  })
  
  datatable_CNA<-reactive({
    inputfile_CNA<-input$mytable_CNA
    if(is.null(inputfile_CNA)){return(NULL)}
    CNA_dtable<-read.csv(inputfile_CNA$datapath, header = input$header, sep = input$sep, stringsAsFactors = F)[,-1]    
    return(CNA_dtable)
  })
  
  datatable_mutation<-reactive({
    inputfile_mutation<-input$mytable_mutation
    if(is.null(inputfile_mutation)){return(NULL)}
    mutation_dtable<-read.csv(inputfile_mutation$datapath, header = input$header, sep = input$sep, stringsAsFactors = F)[,c(-1,-2)]
    return(mutation_dtable)
  })
  
  output$mRNAcontents<-renderTable({
    if (is.null(input$mytable_mrna)){return(NULL)}
    return(datatable_mrna())
  })
  
  output$miRNAcontents <- renderTable({
    if(is.null(input$mytable_mirna)){return(NULL)}
    return(datatable_mirna())
  })
  
  output$CNAcontents <- renderTable({
    if(is.null(input$mytable_CNA)){return(NULL)}
    return(datatable_CNA())
  })
  
  output$Mutation_contents <- renderTable({
    if(is.null(input$mytable_mutation)){return(NULL)}
    return(datatable_mutation())
  })
  
  output$mRNAsummary<-renderPrint({
    if (is.null(input$mytable_mrna)){return("File upload error: No file uploaded")}
    #final_data_gene_expression = calculate_fold_change(datatable_mrna(), "mRNA")
    #return(final_data_gene_expression[1:5, 1:5])
    return(summary(datatable_mrna()))
  })
  
  output$miRNAsummary <- renderPrint({
    if(is.null(input$mytable_mirna)){return("File upload error: No file uploaded")}
    return(summary(datatable_mirna()))
  })
  
  output$CNAsummary <- renderPrint({
    if(is.null(input$mytable_CNA)){return("File upload error: No file uploaded")}
    return(summary(datatable_CNA()))
  })
  
  output$mRNAplot <- renderPlot({    
    mRNA <- datatable_mrna()
    rowsize = 183
    colsize = 121
    mRNA  = mRNA[-1,-1]
    rownames(mRNA) = 1:121
    dataMatrix2 = as.matrix(mRNA)
    DATAMATRIX <- matrix(dataMatrix2, ncol = ncol(dataMatrix2), dimnames = NULL)
    DATAMATRIX = matrix(as.numeric(DATAMATRIX), nrow = nrow(DATAMATRIX), ncol = ncol(DATAMATRIX))
    dataMatrix = DATAMATRIX
    hh <- hclust(dist(dataMatrix))
    dataMatrixOrdered <- dataMatrix[hh$order, ]
    par(mfrow = c(1, 3))
    image(t(dataMatrixOrdered)[, nrow(dataMatrixOrdered):1], xlab = "Patient ID", ylab = "Genes")
    plot(rowMeans(dataMatrixOrdered), colsize:1, , xlab = "Row wise mean expression ", ylab = "Patient ID", pch = 19)
    plot(colMeans(dataMatrixOrdered), xlab = "Genes", ylab = "Column wise mean expression", pch = 19)
    #heatmap.2(dataMatrix, scale = 'row')   
  }, height = 600)
  
  output$miRNAplot <- renderPlot({
    miRNA <- datatable_mirna()
    rowsize = 158
    colsize = 121
    miRNA  = miRNA[-1,-1]
    rownames(miRNA) = 1:121
    dataMatrix2 = as.matrix(miRNA)
    DATAMATRIX <- matrix(dataMatrix2, ncol = ncol(dataMatrix2), dimnames = NULL)
    DATAMATRIX = matrix(as.numeric(DATAMATRIX), nrow = nrow(DATAMATRIX), ncol = ncol(DATAMATRIX))
    dataMatrix = DATAMATRIX
    hh <- hclust(dist(dataMatrix))
    dataMatrixOrdered <- dataMatrix[hh$order, ]
    par(mfrow = c(1, 3))
    image(t(dataMatrixOrdered)[, nrow(dataMatrixOrdered):1], xlab = "Patient ID", ylab = "Genes")
    plot(rowMeans(dataMatrixOrdered), colsize:1, , xlab = "Row wise mean expression ", ylab = "Patient ID", pch = 19)
    plot(colMeans(dataMatrixOrdered), xlab = "Genes", ylab = "Column wise mean expression", pch = 19)
  }, height = 600)
  ############ Score and Analysis #############
  simsScore <- reactive({
    node_names = c("Her_pathway", "CDK4_6", "PLK_AURKA_Kinesins", "ANGIOGENESIS", "ANGIOPOIETINS", "IMMUNO-Modulator", "PI3K", "MET",
                   "MEK", "ERK", "Antiapoptosis", "FGF", "mTOR_AKT_PTEN", "RAS", "TELOMERASE", "IGF_Warburg", "WNT", "PARP", "HDAC",
                   "JAK_STAT", "HEDGEHOG", "NOTCH", "DNA_REPAIR","OTHERS")
    final_score_mrna = calculate_fold_change(datatable_mrna(), "mRNA")
    mrna_score = final_score_mrna[[1]]
    avg_mrna = final_score_mrna[[2]]    
    
    avg_mirna = calculate_fold_change(datatable_mirna(), "miRNA")
        
    matched_expression_mRNA_miRNA = avg_mrna/avg_mirna
    #### Since miRNA NOTCH clumn is NA replace that coulumn with mRNA notch
    matched_expression_mRNA_miRNA[,22] = avg_mrna[, 22]
    mirna_score <- decile_Calculation(matched_expression_mRNA_miRNA, datatable_mrna())
        
    stupid <<- datatable_mrna()
    final_score_CNA = calculate_CNA(datatable_CNA(), datatable_mrna())
    final_score_numeric <<- matrix(NA, nrow = 121, ncol = 24)
    for(j in 1:24){
      final_newscore_CNA = sapply(final_score_CNA[,j], as.numeric)
      final_score_numeric[,j] = final_newscore_CNA
    }
    
    mutation_original = datatable_mutation()
    m_score = calculate_mutation_score(mutation_original, datatable_mrna())
    # Calculate mean #
    mrna_score1 = as.matrix(mrna_score)
    mrna_score2 = matrix(mrna_score1, ncol = 24, dimnames = NULL)
    total_score = mrna_score2 + mirna_score + final_score_numeric + m_score
    total_score = total_score/4
    colnames(total_score) = node_names
    #total_score <- gvisTable(total_score)#, options = list(page = 'enable', height = 300, width = 800))
    #total_score = final_score_CNA
    return(total_score)    
  })


  output$calculateScore <- renderTable({
    return(simsScore())
  })

  
  output$scoreAnalysis <- renderPlot({    
    sims_score <- simsScore()
    histology = histology()
    #print(hist)
    histology <- as.data.frame(histology)
    sims_score <- cbind(histology, sims_score)
    par(mfrow = c(4, 1))
    hc <- hclust(dist(sims_score[, -1], method = "euclidean"), method="ave")
    plot(hc, hang = -1, main = "Score dendrogram for NSCLC cohort", xlab = "Node names") 
    ############ Extract patients with similar histology ########################
    score_table_AC = sims_score[which(sims_score$histology == 'AC'),]
    hc <- hclust(dist(score_table_AC[,-1], method = "euclidean"), method="ave")
    plot(hc, hang = -1, main = "Score dendrogram for NSCLC-AC cohort", xlab = "Node names")
    
    score_table_LCC = sims_score[,-1][which(sims_score$histology == 'LCC'),]
    hc <- hclust(dist(score_table_LCC[, -1], method = "euclidean"), method="ave")
    plot(hc, hang = -1, main = "Score dendrogram for NSCLC-LCC cohort", xlab = "Node names")
    
    score_table_SCC = sims_score[,-1][which(sims_score$histology == 'SCC'),]
    hc <- hclust(dist(score_table_SCC[, -1], method = "euclidean"), method= "ave")
    plot(hc, hang = -1, main = "Score dendrogram for NSCLC-SCC cohort", xlab = "Node names")    
  }, height = 2000)  
  
  ntext <- eventReactive(input$goButton, {
    input$n
  })
  
  output$nText <- renderText({
    ntext()
  })
  
})
