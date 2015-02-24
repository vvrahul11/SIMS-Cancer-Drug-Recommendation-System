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




#decile calculation  # fold_change.csv file was created from final_data_geneexpression file..you can use the original 
## file also but you have to remove the first two columns before startig to use it
data_fold =read.csv("fold_change.csv", as.is = T)
final_data1 = decile_Calculation(data_fold)
write.csv(final_data1, file = "~/Dropbox/winther/Data/score_table_gene_expression.csv")

#decile calculation
data_fold =read.csv("matched_expression_mRNA_miRNA.csv", as.is = T)
final_data2 = decile_Calculation(data_fold)
write.csv(final_data2, file = "~/Dropbox/winther/Data/score_table_mRNA_miRNA.csv")



################# Mean score ###################
mRNA_score = read.csv("score_table_gene_expression.csv", as.is = T)
gene_names = mRNA_score[]
mRNA_score = mRNA_score[,-1:-2]

mRNA_miRNA_score = read.csv("score_table_mRNA_miRNA.csv", as.is = T)
mRNA_miRNA_score = mRNA_miRNA_score[,-1:-2]


# Calculate mean #
total_score = (mRNA_score + mRNA_miRNA_score)/2
colnames(total_score) = node_names
write.csv(total_score, file = "~/Dropbox/winther/Data/total_score.csv")


for(i in 1:24){
  plot(c(1:121), sort(log(data_fold[,i])))
}



