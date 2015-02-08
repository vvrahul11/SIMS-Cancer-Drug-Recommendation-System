inputdata_analysis <- function(dataMatrix1, gene_count){
  #par(mar = rep(0.2, 4))
  myData = matrix(NA, nrow = gene_count, ncol = 121)
  #dataMatrix1 <- read.csv('mir_rev_path-31Jan2014.csv', header = FALSE)
  # CHange gene_count t0 181
  #dataMatrix12 <- read.csv('mRNA_pathway-31Jan2014.csv', header = FALSE)
  
  dataMatrix1  = dataMatrix1[-1:-2,-1]
  rownames(dataMatrix1) = 1:121
  dataMatrix2 = as.matrix(dataMatrix1)
  DATAMATRIX <- matrix(dataMatrix2, ncol = ncol(dataMatrix2), dimnames = NULL)
  DATAMATRIX = matrix(as.numeric(DATAMATRIX), nrow = nrow(DATAMATRIX), ncol = ncol(DATAMATRIX))
  
  image(1:gene_count, 1:121, t(DATAMATRIX)[, nrow(DATAMATRIX):1], xlab = "Genes", ylab = "Patient ID")
  
  dataMatrix = DATAMATRIX
#   par(mar = rep(0.2, 4))
#   image(1:gene_count, 1:121, t(dataMatrix)[, nrow(dataMatrix):1])
#   
  #par(mar = rep(0.2, 4))
  heatmap(dataMatrix)
  
  hh <- hclust(dist(dataMatrix))
  dataMatrixOrdered <- dataMatrix[hh$order, ]
  par(mfrow = c(1, 3))
  image(t(dataMatrixOrdered)[, nrow(dataMatrixOrdered):1], xlab = "Genes", ylab = "Patient ID")
  plot(rowMeans(dataMatrixOrdered), 121:1, , xlab = "Row Mean", ylab = "Row", pch = 15)
  plot(colMeans(dataMatrixOrdered), xlab = "Column", ylab = "Column Mean", pch = 15)
  hist(rowMeans(dataMatrixOrdered), 121:1, , ylab = "Row Mean", xlab = "Row", pch = 15)
  hist(colMeans(dataMatrixOrdered), ylab = "Column", xlab = "Column Mean", pch = 15)

  
  svd1 <- svd(scale(dataMatrixOrdered))
  par(mfrow = c(1, 3))
  image(t(dataMatrixOrdered)[, nrow(dataMatrixOrdered):1])
  plot(svd1$u[, 1], 121:1, , xlab = "Row", ylab = "First left singular vector", 
       pch = 19)
  plot(svd1$v[, 1], xlab = "Column", ylab = "First right singular vector", pch = 19)
  hist(svd1$v[, 1], xlab = "Column", ylab = "First right singular vector", pch = 19)

  
  par(mfrow = c(1, 2))
  plot(svd1$d, xlab = "Column", ylab = "Singular value", pch = 19)
  plot(svd1$d^2/sum(svd1$d^2), xlab = "Column", ylab = "Prop. of variance explained", 
       pch = 19)
  
  
  
  svd1 <- svd(scale(dataMatrixOrdered))
  pca1 <- prcomp(dataMatrixOrdered, scale = TRUE)
  plot(pca1$rotation[, 1], svd1$v[, 1], pch = 19, xlab = "Principal Component 1", 
       ylab = "Right Singular Vector 1")
  abline(c(0, 1))
}

gene_count_mirna = 158
gene_count_mrna = 181
filename_mrna = read.csv('mRNA_pathway-31Jan2014.csv', header = FALSE)
filename_mirna = dataMatrix1 <- read.csv('mir_rev_path-31Jan2014.csv', header = FALSE)
#### Run one function after the other #
inputdata_analysis(filename_mrna, gene_count_mrna)
inputdata_analysis(filename_mirna, gene_count_mirna)
