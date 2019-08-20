#' Calculate kendall values and test
#'
#' f_kendall calculate kendall values, kendallTest test kendall values
#'
#' @param x x
#' @param k k
#' @param geno geno types data
#' @param pheno pheno types data
#' @return p and p values
#' @export
#' @import Kendall

f_kendall <-function(x, k) {

#    if(!require(Kendall)) install.packages("Kendall")
#    library(Kendall)

    ken<-Kendall(x,k) ############### here is changed
    p<-ken$sl

    return(p)
}

kendallTest<-function(geno,pheno) {

        reject<-0
	pvalue<-rep(NA,nrow(geno))

	pvalue<-apply(geno,1,f_kendall,pheno)

	for(n in 1:nrow(geno))
		if(pvalue[n]<0.05)	reject<-reject+1


	print(c("reject count in kendall.R is : ",reject))
	return(pvalue)
}
