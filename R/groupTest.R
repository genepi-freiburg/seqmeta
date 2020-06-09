library(rbgen)
library(seqMeta)

chr = as.numeric(commandArgs(trailingOnly=T)[1])
bgen_path = "/data/studies/06_UKBB/Exome_50k/03_Merge_GWAS_WES/wes/wes_imp_bgen/ukb_wes_efe-chr%CHR%.bgen"
group_file = "/data/studies/06_UKBB/Exome_50k/06_seqMeta/GROUP/03_rareDamaging/group_file.txt"
phenotype_file = "phenotype.txt"
phenotype_col = "egfr_ckdepi_creat"
covariate_cols = "age_crea_serum, U1,U2,U3"
sv_output_file = paste("output/single_variants_chr", chr, ".txt", sep="")
group_output_file = paste("output/group_tests_chr", chr, ".txt", sep="")

# Option for method in skatOMeta function: if lowest p-value in T1 or SKAT
# test < 1E-9, change this to "saddlepoint". From seqMeta package reference
skat_o_method = "integration"



load_genotype = function(chr, bgen_path, snps) {
        bgen_file = gsub("%CHR%", chr, bgen_path)
        print(paste("Loading BGEN file: ", bgen_file, sep=""))

        data = bgen.load(bgen_file, rsids = snps)
        print(paste("Got data for ", length(data$samples), " samples and ", length(data$variants$rsid), " variants", sep=""))

        sum_dosages = 0 * data$data[,,3] + data$data[,,2] + 2 * data$data[,,1]
	if (length(data$variants$rsid) > 1) {
	        sd2 = t(sum_dosages)
		my_geno = cbind(data$samples, sd2)
		colnames(my_geno)[1] = "individual_id"
	} else {
		my_geno = data.frame(individual_id = data$samples, snp = sum_dosages)
		colnames(my_geno)[2] = snps[1]
	}
	print(paste("Genotype data frame: ", dim(my_geno)[1], "x", dim(my_geno)[2], sep=""))
	my_geno
}

prepare_genotype_phenotype_matrices = function(genotype, phenotype) {
	data = merge(phenotype, genotype, by="individual_id")
	print(paste("Merged genotype and phenotype data: ", dim(data)[1], "x", dim(data)[2], sep=""))

	if (ncol(phenotype) + 1 < ncol(data)) {
		geno_matrix = data[, (ncol(phenotype)+1):ncol(data)]
	} else {
		geno_matrix = as.matrix(data[, (ncol(phenotype)+1):ncol(data)])
		colnames(geno_matrix)[1] = colnames(data)[ncol(phenotype)+1]
	}

	for (i in 1:ncol(geno_matrix)) {
		geno_matrix[, i] = as.numeric(as.character(geno_matrix[, i]))
	}

	pheno_matrix = data[, 2:ncol(phenotype)]

        print(paste("Genotype matrix: ", dim(geno_matrix)[1], "x", dim(geno_matrix)[2], 
		"; phenotype matrix: ", dim(pheno_matrix)[1], "x", dim(pheno_matrix)[2], sep=""))
	#print(head(geno_matrix))
	#print(head(pheno_matrix))
	
	list(genotype_matrix = geno_matrix, phenotype_matrix = pheno_matrix)
}

build_snpinfo = function(gene, snps) {
	snp_info = data.frame(
		gene = as.character(rep(gene, length(snps))),
		Name = as.character(snps))
	#print(head(snp_info))
	snp_info
}

calculate_null_model_residuals = function(phenotype_matrix, model_formula) {
	null_model = stats::glm(formula = model_formula, data = phenotype_matrix)
	#print(summary(null_model))
	residuals = stats::residuals(null_model, type = "response")
	print(paste("Calculated ", length(residuals), " residuals for null model.", sep=""))
	residuals
}

process_gene = function(chr, bgen_path, gene, snps, phenotype, model_formula, sv_output_file, group_output_file, skat_o_method = "integration", write_header = F) {
	print(paste("==== Process gene ", gene, " (", length(snps), " variants, chr", chr, ")", sep=""))

	genotype = load_genotype(chr, bgen_path, snps)
	geno_pheno = prepare_genotype_phenotype_matrices(genotype, phenotype)
	snp_info = build_snpinfo(gene, snps)
	#residuals = calculate_null_model_residuals(geno_pheno$phenotype_matrix, model_formula)

	scores = prepScores2(Z = geno_pheno$genotype_matrix,
		formula = model_formula,
		SNPInfo = snp_info,
		data = geno_pheno$phenotype_matrix)
	results = perform_tests(scores, snp_info, skat_o_method)

	print(sv_output_file)

	write_sv_result(results$single_variant, sv_output_file, write_header)
	write_group_result(results, group_output_file, write_header)
}

write_sv_result = function(sv, sv_output_file, write_header) {
	write.table(sv, sv_output_file, append = T, row.names = F, col.names = write_header, sep = "\t", quote = F)
}

write_group_result = function(result, group_output_file, write_header) {
	line = data.frame(gene = result$skat$gene,
			skat_p = result$skat$p,
			skat_Qmeta = result$skat$Qmeta,
			skat_cmaf = result$skat$cmaf,
			skat_nmiss = result$skat$nmiss,
			skat_nsnps = result$skat$nsnps,
			burden_p = result$burden$p,
			burden_beta = result$burden$beta,
			burden_se = result$burden$se,
			burden_cmafTotal = result$burden$cmafTotal,
			burden_cmafUsed = result$burden$cmafUsed,
			burden_nsnpsTotal = result$burden$nsnpsTotal,
			burden_nsnpsUsed = result$burden$nsnpsUsed,
			burden_nmiss = result$burden$nmiss,
			skat_o_p = result$skat_o$p,
			skat_o_pmin = result$skat_o$pmin,
			skat_o_rho = result$skat_o$rho,
			skat_o_cmaf = result$skat_o$cmaf,
			skat_o_nmiss = result$skat_o$nmiss,
			skat_o_nsnps = result$skat_o$nsnps,
			skat_o_errflag = result$skat_o$errflag)
	write.table(line, group_output_file, append = T, row.names = F, col.names = write_header, sep = "\t", quote = F)
}

prepare_phenotype = function(phenotype_file, phenotype_col, individual_col, covariate_cols) {
	phenotypes = read.table(phenotype_file, h=T, sep="\t")
	print(paste("Got phenotypes for ", nrow(phenotypes), " individuals.", sep=""))

	pheno_na = which(is.na(phenotypes[, phenotype_col]))
	print(paste("Excluding ", length(pheno_na), " rows because of missing values for: ", phenotype_col, sep=""))
	phenotypes = phenotypes[-pheno_na,]

	covariates_arr = trimws(unlist(strsplit(covariate_cols, ",", fixed=T)[[1]]))
	for (covar in covariates_arr) {
		covar_na = which(is.na(phenotypes[, covar]))
		if (length(covar_na) > 0) {
			print(paste("Excluding ", length(covar_na), " additional rows because of missing values for: ", covar, sep=""))
			phenotypes = phenotypes[-covar_na,]
		} else {
			print(paste("No missing values for: ", covar, sep=""))
		}
	}

	phenotypes2 = data.frame(individual_id = phenotypes[, individual_col], 
		phenotype = phenotypes[, phenotype_col])
	colnames(phenotypes2)[2] = phenotype_col
	phenotypes2 = cbind(phenotypes2, phenotypes[, covariates_arr])
	print(paste("Phenotype data frame: ", dim(phenotypes2)[1], "x", dim(phenotypes2)[2], sep=""))
	print(summary(phenotypes2))
	phenotypes2
}

prepare_formula = function(phenotype_col, covariate_cols) {
        covariates_arr = trimws(unlist(strsplit(covariate_cols, ",", fixed=T)[[1]]))
	covars_str = paste(covariates_arr, sep="", collapse=" + ")
	formula_str = paste(phenotype_col, " ~ ", covars_str, sep="")
	print(paste("formula as string: ", formula_str, sep=""))
	as.formula(formula_str)
}

perform_tests = function(scores, snp_info, skat_o_method) {
	single_variant = singlesnpMeta(scores, SNPInfo = snp_info)
	#print(format(head(single_variant), digits = 2))

	print("SKAT test")
	skat = skatMeta(scores, SNPInfo = snp_info)
	print(format(head(skat), digits = 2))

	print("Burden test")
	burden = burdenMeta(scores, SNPInfo = snp_info)
	print(format(head(burden), digits=2))

	print("SKAT-O test")
	#meta.skato.results <- skatOMeta(scores, rho=seq(0,1,length=11),
	skat_o = skatOMeta(scores, SNPInfo=snp_info, burden.wts = function(maf) { dbeta(maf, 1, 25) }, method = skat_o_method)
	print(format(head(skat_o), digits=2))

	list(single_variant = single_variant, burden = burden, skat = skat, skat_o = skat_o)
}

print("Clean previous output files")
file.remove(sv_output_file)
file.remove(group_output_file)

print("======= PREPARING PHENOTYPES")

phenotype = prepare_phenotype(phenotype_file, phenotype_col, "individual_id", covariate_cols)
model_formula = prepare_formula(phenotype_col, covariate_cols)

print("======= PROCESS GROUP FILE")

con = file(group_file, "r")
write_header = T
while (TRUE) {
	line = readLines(con, n = 1)
	if (length(line) == 0) {
		break
	}

	gene_and_snps_arr = unlist(strsplit(line, "\t", fixed=T)[[1]])
	gene = gene_and_snps_arr[1]
	snps = gene_and_snps_arr[2:length(gene_and_snps_arr)]

	# TODO chromosome might not be coded in variant identifier
	line_chr = unlist(strsplit(snps[1], ":", fixed=T)[[1]])[1]
	if (chr == line_chr) {
		process_gene(chr, bgen_path, gene, snps, phenotype, model_formula, sv_output_file, group_output_file, skat_o_method, write_header)
		write_header = F
	}
}
close(con)

print("Finished")
