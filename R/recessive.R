#' Create a recessive data
#'
#' recessive function create a recessive pheno and geno type data
#'
#' @param n number of data
#' @param samp number of samples
#' @param mu average
#' @param sigma sigma
#' @param slope slope
#' @return geno and pheno type data and list
#' @export

recessive <- function(n, samp, mu, sigma, slope)
{
	mutation <-0.45  # default value

        change <-0

	value <- c(0,0,slope)
	if(!is.matrix(value)) value <- as.matrix(value)

	## Allocate matrix of arbitrary values by loci
	geno <- matrix(data=1, nrow=n, ncol=samp)
	start_geno <- matrix(NA, nrow=1, ncol=samp)
	pheno <- matrix(data=1, nrow=samp, ncol=1)


	start_geno[1,] <-sample(1:3,samp, replace=TRUE,prob=c(0.25,0.5,0.25))

	pheno<- value[start_geno[1, ],1]
	pheno <- pheno+rnorm(n=samp, mean=mu, sd=sigma)

	start_geno <- start_geno - 1 ## convert 1:3 to 0:2, i.e., number of A1 alleles

	for( i in 1:n)  ## data copy
        	geno[i,] <- start_geno[1,]

        change <- mutation*samp
	for(i in 1:n) {
	    idx <- sample(1:samp, change, replace=F)
	    for(j in 1:samp) { ifelse(j==idx,ifelse(geno[i,j]==0 | geno[i,j]==2, geno[i,j]<-1, geno[i,j]<-3),geno[i,j]<-geno[i,j])}
	}

        for(i in 1:n) {
	    for(j in 1:samp) {
		if(geno[i,j]==3) {
		    q<-runif(1,0,1)
		    if(q<0.5) { geno[i,j]<-0} else{ geno[i,j]<-2}
		}
		else {
		    geno[i,j]<-geno[i,j]
		}
	    }
	}

	write.table(geno, file="sample_genotype.txt", row.names=F, col.names=F, quote=F,sep="\t")
        write.table(pheno, file="sample_phenotype_res.txt", row.names=F, col.names=F, quote=F,sep="\t")

	return (list(geno=geno, pheno=pheno))
}
