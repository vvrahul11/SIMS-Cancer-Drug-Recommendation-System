library(plyr)

#### Copy number alteration

decile_Calculation1 <- function(data_fold, data_foldChange){
  score_matrix= matrix(NA, 121, 24)
  #Read file of fold change for 121 patients, 24 nodes  
  # Get the patient IDs to a vector
  patientID = data_foldChange[,1][-1]
  
  for( i in 1:24){  
    #sorted_data_fold = sort(data_fold[,i], decreasing = FALSE)
    #decile_calculated<-cut(sorted_data_fold,quantile(sorted_data_fold,(0:10)/10),include.lowest=TRUE)
    X = data_fold[,i]
    decile_calculated<- cut(X,quantile(X,(0:10)/10),include.lowest=TRUE)  
    score_matrix[,i] = decile_calculated
  }
  final_data = cbind(patientID, score_matrix)
  return(final_data)
}

calculate_fold_change1 <- function(data_foldChange, datatype, xxx){
  # Code to calculate the average fold change of nodes
  node_names = c("Her_pathway", "CDK4_6", "PLK_AURKA_Kinesins", "ANGIOGENESIS", "ANGIOPOIETINS", "IMMUNO-Modulator", "PI3K", "MET",
                 "MEK", "ERK", "Antiapoptosis", "FGF", "mTOR_AKT_PTEN", "RAS", "TELOMERASE", "IGF_Warburg", "WNT", "PARP", "HDAC",
                 "JAK_STAT", "HEDGEHOG", "NOTCH", "DNA_REPAIR","OTHERS")
  # Get the patient IDs to a vector
  patientID = data_foldChange[,1][-1]
  # Get node names to a vector
  SIMS=data_foldChange[1,]
  # Remove node names from data_foldChange
  new_data_foldChange=data_foldChange[c(-1),]
  # Create a matrix to store the mean fold change
  final_data= matrix(NA, 121, 24)
  
  # For the 24 nodes and 121 patients, calulate the mean fold change
  for(i in 1:24){
    ### Since Notch node is absent in miRNA put NA in column 22 only for miRNA ###
    if(node_names[i] != 'NOTCH' || datatype != "miRNA"){
      sub1=new_data_foldChange[,which(SIMS== node_names[i])]    
      sub2=matrix(NA,ncol=ncol(sub1),nrow=nrow(sub1))
      for(j in 1:ncol(sub1)){
        sub2[,j]=as.numeric(sub1[,j])
      }
      sub2 = as.matrix(sub2)
      fold_change = rowMeans(sub2)
      final_data[,i] = fold_change      
    }    
  }
  if(datatype == "mRNA"){
    # Add nodes names( as column names) to the mean fold change
    colnames(final_data) = node_names  
    final_score = decile_Calculation1(final_data, xxx)    
    final_score = as.data.frame(final_score)
    #names(final_score) = node_names
    return(list(final_score, final_data))
  }
  
  else{
    return(final_data)    
  }
  
}


calculate_CNA <- function(cna, mrna_full){
  xxx = mrna_full
  mrna = mrna_full[-1,-1]
  #CNA_MRNA_matrix = matrix(NA, nrow = 121, ncol = 181 )
  CNA_MRNA_matrix = data.frame(matrix(NA, nrow = 121, ncol = 181 ))
  for(i in 1:121){
    x = as.numeric(mrna[i,])
    y = cna[i,]
    z = as.numeric(revalue(as.character(y), c(N = 1, A = 1.5, D = 0.5)), warn_missing = F)
    prod = x * z
    print(prod)
    CNA_MRNA_matrix[i, ] = prod   
  } 
  mrna_full[2:122, 2:182] = CNA_MRNA_matrix  
  # Code to calculate the average fold change of nodes
  node_names = c("Her_pathway", "CDK4_6", "PLK_AURKA_Kinesins", "ANGIOGENESIS", "ANGIOPOIETINS", "IMMUNO-Modulator", "PI3K", "MET",
                 "MEK", "ERK", "Antiapoptosis", "FGF", "mTOR_AKT_PTEN", "RAS", "TELOMERASE", "IGF_Warburg", "WNT", "PARP", "HDAC",
                 "JAK_STAT", "HEDGEHOG", "NOTCH", "DNA_REPAIR","OTHERS")
  
  final_data_cnv = calculate_fold_change1(mrna_full, "mRNA", xxx)  
  cna_score = as.matrix(final_data_cnv[[1]][,-1])
  #print(cna_score)
  return(cna_score)
}

calculate_mutation_score <- function(mutation_original, mrna){
  #### Mutation information
  mutation = apply(mutation_original, 2, function(x) gsub("^$|^ $", 0, x))
  mutation = apply(mutation, 2, function(x) gsub("WT", 0, x))
  mutation = apply(mutation, 2, function(x) gsub("wt", 0, x))
  mutation_matrix = matrix(NA, nrow = 123, ncol = 6)
  for(i in 1:dim(mutation)[2]){
    mutation[,i] <- as.character(mutation[,i])
    mutation[,i][mutation[,i] != 0] <- 1  
  }
  ####### Mutation score #########
  nodes = c("Her_pathway","CDK4_6","PLK_AURKA_Kinesins","ANGIOGENESIS","ANGIOPOIETINS","IMMUNO-Modulator","PI3K","MET","MEK","ERK","Antiapoptosis","FGF","mTOR_AKT_PTEN","RAS","TELOMERASE","IGF_Warburg","WNT","PARP","HDAC","JAK_STAT","HEDGEHOG","NOTCH","DNA_REPAIR","OTHERS")
  mutation_score = matrix(0, nrow = 121, ncol = 24)
  colnames(mutation_score) = nodes
  mutated_genes = colnames(mutation)
  mutation_score[,14]  = mutation[,1]
  mutation_score[,1]  = mutation[,2]
  mutation_score[,11]  = mutation[,3]
  mutation_score[,14]  = mutation[,4]
  mutation_score[,1]  = mutation[,5]
  ### this should have been for P53 but since its expression is absent
  ### another gene in the same pathway is selected bcz the score wil finally effect the node 
  mutation_score[,5]  = mutation[,6] 
  mutation_score[mutation_score == 1] <- 10
  mutation_score_numeric = matrix(NA, nrow = 121, ncol = 24)  
  for(i in 1:24){
    column_numeric_mutation = sapply(mutation_score[,i], as.numeric)
    mutation_score_numeric[,i] = column_numeric_mutation
  }
  return(mutation_score_numeric)
}

decile_Calculation <- function(data_fold, data_foldChange){
  score_matrix= matrix(NA, 121, 24)
  #Read file of fold change for 121 patients, 24 nodes
  # Get the patient IDs to a vector
  patientID = data_foldChange[,1][-1]
  
  for( i in 1:24){  
    #sorted_data_fold = sort(data_fold[,i], decreasing = FALSE)
    #decile_calculated<-cut(sorted_data_fold,quantile(sorted_data_fold,(0:10)/10),include.lowest=TRUE)
    X = data_fold[,i]
    decile_calculated<- cut(X,quantile(X,(0:10)/10),include.lowest=TRUE)  
    score_matrix[,i] = decile_calculated
  }
  #print(dim(score_matrix))
  #final_data3 = cbind(patientID, score_matrix)
  return(score_matrix)
}

calculate_fold_change <- function(data_foldChange, datatype){
  
  # Code to calculate the average fold change of nodes
  node_names = c("Her_pathway", "CDK4_6", "PLK_AURKA_Kinesins", "ANGIOGENESIS", "ANGIOPOIETINS", "IMMUNO-Modulator", "PI3K", "MET",
                 "MEK", "ERK", "Antiapoptosis", "FGF", "mTOR_AKT_PTEN", "RAS", "TELOMERASE", "IGF_Warburg", "WNT", "PARP", "HDAC",
                 "JAK_STAT", "HEDGEHOG", "NOTCH", "DNA_REPAIR","OTHERS")
  # Get the patient IDs to a vector
  patientID = data_foldChange[,1][-1]
  # Get node names to a vector
  SIMS.node_names=data_foldChange[1,]
  # Remove node names from data_foldChange
  new_data_foldChange=data_foldChange[c(-1),]
  # Create a matrix to store the mean fold change
  final_data= matrix(NA, nrow(data_foldChange)-1, length(node_names))
  
  # For the 24 nodes and 121 patients, calulate the mean fold change
  for(i in 1:24){
    ### Since Notch node is absent in miRNA put NA in column 22 only for miRNA ###
    if(node_names[i] != 'NOTCH' || datatype != "miRNA"){
      sub1=new_data_foldChange[,which(SIMS.node_names== node_names[i])]    
      sub2=matrix(NA,ncol=ncol(sub1),nrow=nrow(sub1))
      for(j in 1:ncol(sub1)){
        sub2[,j]=as.numeric(sub1[,j])
      }
      sub2 = as.matrix(sub2)
      fold_change = rowMeans(sub2)
      final_data[,i] = fold_change      
    }    
  }
  
  if(datatype == "mRNA"){
    # Add nodes names( as column names) to the mean fold change
    colnames(final_data) = node_names  
    final_score = decile_Calculation(final_data, data_foldChange)
    final_score = as.data.frame(final_score)
    names(final_score) = node_names
    out_object <- list(final_score=final_score, final_data=final_data);
    return(out_object)
  }
  
  else{
    return(final_data)    
  }
  
}

# Define input parameters for main analysis (done in "calc_SIMS")
input = list(
  mrna=list(
    datapath="SIMS.input.mRNA.csv"
  ),
  mirna=list(
    datapath="SIMS.input.miRNA.csv"
  ),
  CNA=list(
    datapath="SIMS.input.CNA.csv"
  ),
  mutation=list(
    datapath="SIMS.input.mutation.csv"
  ),
  node_names = c("Her_pathway", "CDK4_6", "PLK_AURKA_Kinesins", "ANGIOGENESIS", "ANGIOPOIETINS", "IMMUNO-Modulator", "PI3K", "MET",
                 "MEK", "ERK", "Antiapoptosis", "FGF", "mTOR_AKT_PTEN", "RAS", "TELOMERASE", "IGF_Warburg", "WNT", "PARP", "HDAC",
                 "JAK_STAT", "HEDGEHOG", "NOTCH", "DNA_REPAIR","OTHERS")
  
)

calc_SIMS<-function(input) { 
  
  output=list()
  ############ Read data from files #############
  # Generate a summary of the data
  input$mrna_dtable<-read.csv(input$mrna$datapath, stringsAsFactors = F)    # header for header file, sep for separator
  input$Histology <- c("AC","AC","SCC","AC","AC","SCC","AC","AC","SCC","SCC","LCC","AC","SCC","AC","SCC","AC","AC","AC","Other SCLC","LCC","AC","LCC","LCC","AC","AC","SCC","SCC","SCC","AC","SCC","AC","AC","Other ADEC","SCC","SCC","AC","SCC","LCC","AC","AC","SCC","AC","AC","AC","SCC","AC","SCC","SCC","SCC","AC","SCC","AC","SCC","AC","AC","SCC","AC","AC","AC","AC","SCC","SCC","AC","AC","SCC","SCC","AC","LCC","SCC","AC","AC","AC","LCC","SCC","SCC","SCC","SCC","AC","LCC","SCC","AC","SCC","AC","SCC","SCC","SCC","AC","SCC","LCC","SCC","AC","AC","AC","AC","LCC","AC","Other ADEC","AC","SCC","SCC","AC","LCC","SCC","AC","SCC","AC","SCC","LCC","AC","SCC","AC","SCC","AC","SCC","SCC","SCC","SCC","SCC","SCC","AC","AC")
  input$mirna_dtable<-read.csv(input$mirna$datapath, stringsAsFactors = F)    
  input$CNA_dtable<-read.csv(input$CNA$datapath, stringsAsFactors = F)    
  input$mutation_dtable<-read.csv(input$mutation$datapath, stringsAsFactors = F)    
  
  if ("part2" != "store data for later use in output") {
    if (is.null(input$mrna_dtable)) {
      output$mRNAcontents<-NULL
    } else {
      output$mRNAcontents<-input$mrna_dtable
    }
    
    if (is.null(input$mirna_dtable)) {
      output$miRNAcontents<-NULL
    } else {
      output$miRNAcontents<-input$mirna_dtable
    }
    
    if (is.null(input$CNA_dtable)) {
      output$CNAcontents<-NULL
    } else {
      output$CNAcontents<-input$CNA_dtable
    }
    
    if (is.null(input$mutation_dtable)) {
      output$Mutation_contents<-NULL
    } else {
      output$Mutation_contents<-input$mutation_dtable
    }
  }
  
  # output$mRNAplot <- renderPlot({    
  #   mRNA <- datatable_mrna()
  #   rowsize = 183
  #   colsize = 121
  #   mRNA  = mRNA[-1,-1]
  #   rownames(mRNA) = 1:121
  #   dataMatrix2 = as.matrix(mRNA)
  #   DATAMATRIX <- matrix(dataMatrix2, ncol = ncol(dataMatrix2), dimnames = NULL)
  #   DATAMATRIX = matrix(as.numeric(DATAMATRIX), nrow = nrow(DATAMATRIX), ncol = ncol(DATAMATRIX))
  #   dataMatrix = DATAMATRIX
  #   hh <- hclust(dist(dataMatrix))
  #   dataMatrixOrdered <- dataMatrix[hh$order, ]
  #   par(mfrow = c(1, 3))
  #   image(t(dataMatrixOrdered)[, nrow(dataMatrixOrdered):1], xlab = "Patient ID", ylab = "Genes")
  #   plot(rowMeans(dataMatrixOrdered), colsize:1, , xlab = "Row wise mean expression ", ylab = "Patient ID", pch = 19)
  #   plot(colMeans(dataMatrixOrdered), xlab = "Genes", ylab = "Column wise mean expression", pch = 19)
  #   #heatmap.2(dataMatrix, scale = 'row')   
  # }, height = 600)
  # 
  # output$miRNAplot <- renderPlot({
  #   miRNA <- datatable_mirna()
  #   rowsize = 158
  #   colsize = 121
  #   miRNA  = miRNA[-1,-1]
  #   rownames(miRNA) = 1:121
  #   dataMatrix2 = as.matrix(miRNA)
  #   DATAMATRIX <- matrix(dataMatrix2, ncol = ncol(dataMatrix2), dimnames = NULL)
  #   DATAMATRIX = matrix(as.numeric(DATAMATRIX), nrow = nrow(DATAMATRIX), ncol = ncol(DATAMATRIX))
  #   dataMatrix = DATAMATRIX
  #   hh <- hclust(dist(dataMatrix))
  #   dataMatrixOrdered <- dataMatrix[hh$order, ]
  #   par(mfrow = c(1, 3))
  #   image(t(dataMatrixOrdered)[, nrow(dataMatrixOrdered):1], xlab = "Patient ID", ylab = "Genes")
  #   plot(rowMeans(dataMatrixOrdered), colsize:1, , xlab = "Row wise mean expression ", ylab = "Patient ID", pch = 19)
  #   plot(colMeans(dataMatrixOrdered), xlab = "Genes", ylab = "Column wise mean expression", pch = 19)
  # }, height = 600)
  
  ############ Score and Analysis #############
  output$node_names = c("Her_pathway", "CDK4_6", "PLK_AURKA_Kinesins", "ANGIOGENESIS", "ANGIOPOIETINS", "IMMUNO-Modulator", "PI3K", "MET",
                        "MEK", "ERK", "Antiapoptosis", "FGF", "mTOR_AKT_PTEN", "RAS", "TELOMERASE", "IGF_Warburg", "WNT", "PARP", "HDAC",
                        "JAK_STAT", "HEDGEHOG", "NOTCH", "DNA_REPAIR","OTHERS")
  final_score_mrna = calculate_fold_change(input$mrna_dtable, "mRNA")

  print(names(final_score_mrna))
  stop("p 286")  
  output$simsScore <- {
    mrna_score = final_score_mrna[[1]]
    avg_mrna = final_score_mrna[[2]]    
    
    avg_mirna = calculate_fold_change(input$mirna_dtable, "miRNA")
    
    print(dim(avg_mirna))
    stop("tag 298")
    print(names(avg_mirna))
    matched_expression_mRNA_miRNA = avg_mrna/avg_mirna
    
    #### Since miRNA NOTCH clumn is NA replace that coulumn with mRNA notch
    matched_expression_mRNA_miRNA[,22] = avg_mrna[, 22]
    mirna_score <- decile_Calculation(matched_expression_mRNA_miRNA, input$mrna_dtable)
    
    stupid <<- input$mrna_dtable
    final_score_CNA = calculate_CNA(input$CNA_dtable, input$mrna_dtable)
    final_score_numeric <<- matrix(NA, nrow = 121, ncol = 24)
    for(j in 1:24){
      print(i)
      final_newscore_CNA = sapply(final_score_CNA[,j], as.numeric)
      final_score_numeric[,j] = final_newscore_CNA
    }
    
    mutation_original = intput$mutation_dtable
    m_score = calculate_mutation_score(mutation_original, input$mrna_dtable)
    # Calculate mean #
    mrna_score1 = as.matrix(mrna_score)
    mrna_score2 = matrix(mrna_score1, ncol = 24, dimnames = NULL)
    total_score = mrna_score2 + mirna_score + final_score_numeric + m_score
    total_score = total_score/4
    colnames(total_score) = node_names
    #total_score <- gvisTable(total_score)#, options = list(page = 'enable', height = 300, width = 800))
    #total_score = final_score_CNA
    return(total_score)    
  }
  
  
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
  
}

output <- calc_SIMS(input)