#' adjgeno
#'
#' adjgeno function is ~
#'
#' @param geno geno
#' @param evec evec
#' @param pca_num pca number
#' @return geno list?
#' @export

# evec <- c()

adjgeno<-function(geno,evec,pca_num) {
	evec<-as.matrix(evec)

	geno1<-geno%*%evec[,1:pca_num]%*%t(evec[,1:pca_num])

	geno <- geno-geno1;

	return(geno)
}
