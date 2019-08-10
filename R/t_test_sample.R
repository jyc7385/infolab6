#' Calculate t values and test
#'
#' f_kendall calculate t values, kendallTest test t values
#'
#' @param x x
#' @param k k
#' @param geno geno types data
#' @param pheno pheno types data
#' @return p and p values
#' @export

f_t_test <- function(x,k) {
	m1 <-lm(k~x)
	p<-summary(m1)$coefficients[2,4]

	return(p)
}

t_test<-function(geno,pheno) {

	reject<-0
	pvalue<-rep(NA,nrow(geno))


	pvalue<-apply(geno,1,f_t_test,pheno)

	for(n in 1:nrow(geno))
		if(pvalue[n]<0.05)	reject<-reject+1

	print(c("reject count in t-test is : ",reject))
	return(pvalue)
}
