#### Copy number alteration

decile_Calculation <- function(data_fold){
  score_matrix= matrix(NA, 121, 24)
  #Read file of fold change for 121 patients, 24 nodes
  data_foldChange =read.csv("mRNA_pathway-31Jan2014.csv", as.is = T)
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


calculate_fold_change <- function(data_foldChange, datatype){
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
    final_score = decile_Calculation(final_data)
    final_score = as.data.frame(final_score)
    #names(final_score) = node_names
    return(list(final_score, final_data))
  }
  
  else{
    return(final_data)    
  }
  
}


calculate_CNA <- function(cna, mrna_full){
  mrna = mrna_full[-1,-1]
  #CNA_MRNA_matrix = matrix(NA, nrow = 121, ncol = 181 )
  CNA_MRNA_matrix = data.frame(matrix(NA, nrow = 121, ncol = 181 ))
  for(i in 1:121){
    x = as.numeric(mrna[i,])
    y = cna[i,]
    z = as.numeric(revalue(as.character(y), c(N = 1, A = 1.5, D = 0.5)), warn_missing = F)
    prod = x * z
    CNA_MRNA_matrix[i, ] = prod   
  } 
  mrna_full[2:122, 2:182] = CNA_MRNA_matrix
  £write.csv(CNA_MRNA_matrix, "copy_number_expression.csv", row.names = F, quote = F)
  # Code to calculate the average fold change of nodes
  node_names = c("Her_pathway", "CDK4_6", "PLK_AURKA_Kinesins", "ANGIOGENESIS", "ANGIOPOIETINS", "IMMUNO-Modulator", "PI3K", "MET",
                 "MEK", "ERK", "Antiapoptosis", "FGF", "mTOR_AKT_PTEN", "RAS", "TELOMERASE", "IGF_Warburg", "WNT", "PARP", "HDAC",
                 "JAK_STAT", "HEDGEHOG", "NOTCH", "DNA_REPAIR","OTHERS")
  #data_cnv =read.csv("copy_number_expression.csv", as.is = T)
  final_data_cnv = calculate_fold_change(mrna_full, "mRNA")
  #write.csv(final_data_cnv[,-1], "fold_change.csv", row.names = F, quote = F)
  #data_fold =read.csv("fold_change.csv", as.is = T)
  #final_data1 = decile_Calculation(data_fold)
  cna_score = as.matrix(final_data_cnv[[1]][,-1])
  print(cna_score)
  return(cna_score)
}
