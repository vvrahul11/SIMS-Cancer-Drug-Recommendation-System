library(gplots)
score_analysis <- function(DATAMATRIX, rowsize, colsize){
  dataMatrix = DATAMATRIX
  par(mar = rep(0.2, 4))
  image(1:rowsize, 1:colsize, t(dataMatrix)[, nrow(dataMatrix):1], xlab = "Patient ID", ylab = 'Genes')
  
  par(mar = rep(0.2, 4))
  heatmap.2(dataMatrix, scale = 'row')
  
  hh <- hclust(dist(dataMatrix))
  dataMatrixOrdered <- dataMatrix[hh$order, ]
  par(mfrow = c(1, 3))
  image(t(dataMatrixOrdered)[, nrow(dataMatrixOrdered):1], xlab = "Patient ID", ylab = "Genes")
  plot(rowMeans(dataMatrixOrdered), colsize:1, , xlab = "Row wise mean expression ", ylab = "Patient ID", pch = 19)
  plot(colMeans(dataMatrixOrdered), xlab = "Column", ylab = "Column wise mean expression", pch = 19)
  
  
  svd1 <- svd(scale(dataMatrixOrdered))
  par(mfrow = c(1, 3))
  image(t(dataMatrixOrdered)[, nrow(dataMatrixOrdered):1])
  plot(svd1$u[, 1], colsize:1, , xlab = "Row", ylab = "First left singular vector", 
       pch = 19)
  plot(svd1$v[, 1], xlab = "Column", ylab = "First right singular vector", pch = 19)
  
  
  par(mfrow = c(1, 2))
  plot(svd1$d, xlab = "Column", ylab = "Singular value", pch = 19)
  plot(svd1$d^2/sum(svd1$d^2), xlab = "Column", ylab = "Prop. of variance explained", 
       pch = 19)
  
  
  
  svd1 <- svd(scale(dataMatrixOrdered))
  pca1 <- prcomp(dataMatrixOrdered, scale = TRUE)
  plot(pca1$rotation[, 1], svd1$v[, 1], pch = 19, xlab = "Principal Component 1", 
       ylab = "Right Singular Vector 1")
  abline(c(0, 1))
  plot(pca1$rotation[,1], pca1$rotation[,2], xlab = "Principal Component 1",
       ylab = "Principal Component 3")
  
}


