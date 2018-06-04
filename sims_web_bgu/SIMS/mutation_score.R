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
