library(  quadprog)
library(Matrix)

dataX <- matrix(c(1,1,2,2,3,3,1,2,1,3,2,3),6,2)
dataY <- c(rep(-1,3),rep(1,3))


findWsBs <- function(dataX,dataY) {
	m <- length(dataY)
	dvec <- rep(1,m)
	Dmat <- matrix(0,m,m)
	for (i in 1:m) 
		for (j in 1:m) 
			Dmat[i,j] = (dataY[i]*dataY[j])*(dataX[i,]%*%dataX[j,])
	Amat <- t ( rbind(dataY,diag(m)) )
	bvec <- rep(0,m+1)
	Dmat <- nearPD(Dmat)$mat
	sol <- solve.QP(Dmat, dvec, Amat, meq=1)
	alpha <- sol$solution
	w <- colSums((alpha*dataY)*dataX)
	w <- w / sqrt(sum(w^2))
	dataXPos <- dataX[which(dataY==1),]
	dataXNeg <- dataX[which(dataY==-1),]
	b = - (  min(dataXPos %*% w) + max(dataXNeg %*% w) ) / 2
	return (list(w=w, b=b))
}


ans <- findWsBs (dataX,dataY)
w <- ans$w
b <- ans$b 

plot(dataX[,1],dataX[,2],col=(dataY+2),pch=20)
x1 <- seq(min(dataX[,1]), max(dataX[,1]), length.out=100)
x2 <- - ( (w[1]/w[2]) * x1 ) - (b/w[2])
lines(x1,x2)

testDataX <- matrix(c(0,4,0,4),2,2)

SVM <- function(dataX, dataY, testDataX) {
	ans <- findWsBs (dataX,dataY)
	w <- ans$w
	b <- ans$b
	return ( sign( (testDataX %*% w) + b ) )
}

pred <- SVM(dataX, dataY, testDataX)
