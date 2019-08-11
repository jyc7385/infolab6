#' Use Data and Calculate
#'
#' GWAS_main function use data and calculate values
#'
#' @param Y phenotypes data
#' @param G genotypes data
#' @param pc_num the number of PCA
#' @return many data
#' @export


# Input : phenotypeData, genotypeData(vcf or hapmap), the number of PCA

if(!require(Kendall)) install.packages("Kendall")
library(Kendall)

GWAS_main<-function(Y,G,pc_num) {

	if(is.null(Y))	stop("Phenotypes must exist.")
	if(is.null(G))	stop("Genotypes must exist.")

	# genotype data loading
	geno_num <-as.matrix(G)
	geno_pos <-1:nrow(geno_num)
	geno_chr <-rep(1,nrow(geno_num))
	bp <- 1:nrow(geno_num)

	# phenotype data loading
	pheno_num<-as.matrix(Y)
	print("Genotype encoding complete.")

	norm_geno <- scale(G, center=T, scale=T) # genotype matirx normalization
	print("Genotype normalization complete.")

	#PCA
	pc <-prcomp(na.omit(norm_geno))

	print("eigenvalue sum is ")
	print(sum(pc$sdev*pc$sdev))

	pca_eval <-pc$sdev #eigen values
	write.table(pca_eval,file="PCA_eigenvalues.txt",row.names=F,col.names=F,quote=F,append=F,sep="\t")
	pca_evec<-pc$rotation[,1:pc_num] #eigen vector
	write.table(pca_evec,file="PCA_eigenvectors.txt",row.names=F,col.names=F,quote=F,append=F,sep="\t")
	print("PCA complete.")

	adj_geno<-adjgeno(geno=geno_num,evec=pca_evec,pc_num)
        write.table(adj_geno,file="adjusted_genotype.txt",row.names=F,col.names=F,quote=F,append=F,sep="\t")

	print("Adjusted Genotype matrix complete.")

	kpv<-kendallTest(geno=adj_geno,pheno=pheno_num)
	ken_pv <-cbind(geno_chr,geno_pos,kpv, bp) # chr snp  p + bp
        write.table(ken_pv,file="kendall_test.txt",row.names=F,col.names=F,quote=F,append=F)
	print("Kendall Test complete.")

	tpv<- t_test(geno=adj_geno,pheno=pheno_num)
	t_pv <-data.frame(geno_chr,geno_pos, tpv)
        write.table(t_pv,file="t_test.txt",row.names=F,col.names=F,quote=F, append=F)
	print("t-test Test complete.")

	adj_kpv<-p.adjust(kpv,method="BH")
        print(c("reject count using adjusted p-value is ",sum(adj_kpv<0.05)))

	adj_tpv<-p.adjust(tpv,method="BH")
        print(c("reject count using adjusted p-value is ",sum(adj_tpv<0.05)))

	adj_kpv <-data.frame(geno_chr,geno_pos, adj_kpv)
	adj_tpv <-data.frame(geno_chr,geno_pos, adj_tpv)

	write.table(adj_kpv,file="adj_kendall_test.txt",sep="\t",row.names=F,col.names=F,quote=F,append=F)
	write.table(adj_tpv,file="adj_t_test.txt",sep="\t",row.names=F,col.names=F,quote=F,append=F)
	print("FDR Test complete.")
}
