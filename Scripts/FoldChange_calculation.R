decile_Calculation <- function(data_fold){
  final_data = final_data[, -1:-2]
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
  print(class(data_foldChange))
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
  # Add nodes names( as column names) to the mean fold change
  colnames(final_data) = node_names
  # Add patientIDs to the calculated mean fold change
  final_data = cbind(patientID, final_data)
  final_data1 = decile_Calculation(final_data)
  return(final_data)
}
  


