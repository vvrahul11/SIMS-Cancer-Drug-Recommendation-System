
####### Create Dendrograms for the score table #########
score_table = read.csv("total_score.csv")
score_table_new = score_table[,-1]
#score_table_new = score_table[,-1]
hc <- hclust(dist(score_table_new[1:121, 1:24], method = "euclidean"), method="ave")
plot(hc, hang = -1)

############ Patient histology count ##########
# AC = 56
# LCC = 12
# SCC = 50
# Others = 3   # ADEC = 33, 97  # SCLC = 19

############ Extract patients with similar histology ########################
score_table_new_AC = score_table[which(score_table$Histology == 'AC'),]
hc <- hclust(dist(score_table_new_AC[,-1], method = "euclidean"), method="ave")
plot(hc, hang = -1)

score_table_new_LCC = score_table[which(score_table$Histology == 'LCC'),]
hc <- hclust(dist(score_table_new_LCC[, -1], method = "euclidean"), method="ave")
plot(hc, hang = -1)

score_table_new_SCC = score_table[which(score_table$Histology == 'SCC'),]
hc <- hclust(dist(score_table_new_SCC[, -1], method = "euclidean"), method="ave")
plot(hc, hang = -1)


########## Standing out unique dendrograms ########
# 
# score_table[61,][1] # Histologies are SCC AC AC
# score_table[64,][1]
# score_table[43,][1]
# 
# score_table[94,][1] # All 3s histology is AC
# score_table[16,][1]
# score_table[81,][1]   
# 
# score_table[76,][1] # All 3s histology is AC
# score_table[56,][1]
# score_table[79,][1]   
# score_table[102,][1] 
