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
