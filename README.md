if(!require(devtools)) install.packages("devtools")

devtools::install_github("jyc7385/infolab6")

library(infolab6)

add_data <- additive(10000, 100, 0, 1, 0.2)
GWAS_main(Y = add_data$pheno, G = add_data$geno, pc_num = 5)

# dom_data <- dominant(10000, 100, 0, 1, 0.2)
# GWAS_main(Y = dom_data$pheno, G = dom_data$geno, pc_num = 5)

# rec_data <- recessive(10000, 100, 0, 1, 0.2)
# GWAS_main(Y = rec_data$pheno, G = rec_data$geno, pc_num = 5)

bio_shiny()
