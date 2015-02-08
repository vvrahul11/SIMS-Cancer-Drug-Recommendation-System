# Read mRNA and miRNA fold change values
data_foldChange1 =read.csv("mRNA_pathway-31Jan2014.csv", as.is = T)
data_foldChange2 =read.csv("mir_rev_path-31Jan2014-2.csv", as.is = T)

# get the gene names and remove the first patient IDs from the column
gene_names1 = names(data_foldChange1)[-1] 
gene_names2 = names(data_foldChange2)[-1]

# Number of genes in mRNA and miRNA fold change
length(gene_names1)  # 181 
length(gene_names2)  # 158

# Number of genes similar in both is 155
genes_index = which(gene_names1 %in% gene_names2)
genes_similar_in_both = gene_names1[which(gene_names1 %in% gene_names2)]
#print(genes_similar_in_both)
print(length(genes_similar_in_both))

# Genes present in miRNA but not in mRNA - position 47 121 148 - names - FGF9 PARP1 HSP90AA1
which(!(gene_names2 %in% gene_names1))
#gene_names2[which(!(gene_names2 %in% gene_names1))]

# Genes present in mRNA but not in miRNA - position 3  36  40  41  48  67  82  85  93 101 121 122 123 124 125 126 127 143 147 148 151 153 163 169 174 179
#  names - "FLT4"   "MAPK1"  "FGF10"  "FGF11"  "FGF5"   "SUFU"   "IGF1R"  "INSR"   "JAK1"   "MAP2K3" "ADAM17" "APH1A"  "JAG1"  
# "NCSTN"  "NOTCH1" "PSEN1"  "SRRT"   "PIK3CB" "PIK3R2" "PIK3R3" "BORA"   "KIF11"  "CCND1"  "PTGES3" "CTNNB1" "WNT1" 
which(!(gene_names1 %in% gene_names2))
# gene_names1[which(!(gene_names1 %in% gene_names2))] 

###### Plot the count of genes and drugs in nodes ##########
genes_nodes  = c(13, 10, 5, 12, 8, 5, 10, 4, 8, 4, 5, 18, 10, 5, 6, 7, 10, 10, 5, 6, 7, 7, 7, 4)
node_names = c("Her", "CDK", "PLK", "A-GEN", "A-POI", "IMMUN", "PI3K", "MET",
               "MEK", "ERK", "Antiapo", "FGF", "mTOR", "RAS", "TELO", "IGF", "WNT", "PARP", "HDAC",
               "JAK", "HEDGE", "FGF", "DNA_RE","OTH")
drugs = c(2,1,1,2,0,7,3,5,2,0,2,3,7,4,0,2,1,3,1,2,1,1,1,4)
par(mfrow = c(2, 1))
barplot(genes_nodes, xlab = "Gene Names", ylab = "Count", col = "blue", main = "Gene count in nodes(No of Genes = 186, No of Nodes = 24)", names.arg = node_names, ylim = c(0, 20))
barplot(drugs, xlab = "Gene Names", ylab = "Count", col = "red", main = "Drug count in nodes(No of Drugs = 58, No of Nodes = 24)", names.arg = node_names, ylim = c(0,10))
