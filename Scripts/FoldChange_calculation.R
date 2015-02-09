calculate_fold_change <- function(node_names, data_foldChange, datatype){
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
  return(final_data)
}
  

#################################################################
######################  Main  ###################################
# Code to calculate the average fold change of nodes
node_names = c("Her_pathway", "CDK4_6", "PLK_AURKA_Kinesins", "ANGIOGENESIS", "ANGIOPOIETINS", "IMMUNO-Modulator", "PI3K", "MET",
               "MEK", "ERK", "Antiapoptosis", "FGF", "mTOR_AKT_PTEN", "RAS", "TELOMERASE", "IGF_Warburg", "WNT", "PARP", "HDAC",
               "JAK_STAT", "HEDGEHOG", "NOTCH", "DNA_REPAIR","OTHERS")

#Extra nodes = "Resist_CDK4_6_inhibit", "Sensitivity_to_CDK4_6_inhibitors"?

#Read file of fold change for 121 patients, 24 nodes
data_foldChange1 =read.csv("mRNA_pathway-31Jan2014.csv", as.is = T)
final_data_gene_expression = calculate_fold_change (node_names, data_foldChange1, "mRNA")
write.csv(final_data_gene_expression, file = "~/Dropbox/winther/Data/final_data_gene_expression.csv")  

data_foldChange2 =read.csv("mir_rev_path-31Jan2014-2.csv", as.is = T)
final_data_miRNA = calculate_fold_change (node_names, data_foldChange2, "miRNA")
write.csv(final_data_miRNA, file = "~/Dropbox/winther/Data/final_data_miRNA.csv")  

mRNA = read.csv("final_data_gene_expression.csv", as.is = T)
mRNA = mRNA[,-1]
mRNA = mRNA[,-1]

miRNA = read.csv("final_data_miRNA.csv", as.is = T)
miRNA = miRNA[,-1]
miRNA = miRNA[,-1]

matched_expression_mRNA_miRNA = mRNA/miRNA
#### Since miRNA NOTCH clumn is NA replace that coulumn with mRNA notch
matched_expression_mRNA_miRNA[,22] = mRNA[, 22]
write.csv(matched_expression_mRNA_miRNA, file = "~/Dropbox/winther/Data/matched_expression_mRNA_miRNA.csv")  


