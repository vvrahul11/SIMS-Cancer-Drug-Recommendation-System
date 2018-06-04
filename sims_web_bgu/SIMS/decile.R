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
