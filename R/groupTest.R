#!/usr/local/R/R-3.6.3/bin/Rscript
print("Loading packages")
suppressPackageStartupMessages(library(rbgen))
suppressPackageStartupMessages(library(seqMeta))
suppressPackageStartupMessages(library(optparse))

parse_options = function() {
  option_list = list( 
    make_option("--chr", help="Chromosome number (for patterns)"),
    make_option("--bgen_path", help="BGEN/BGI file name, %CHR% will be substituted by 'chr'"),
    make_option("--group_file", help="Group file name"),
    make_option("--phenotype_file", help="Phenotype file name"),
    make_option("--phenotype_col", help="Column name of phenotype"),
    make_option("--phenotype_type", help="Phenotype type ('binary'/'quantitative'), default: quantitative", default="quantitative"),
    make_option("--covariate_cols", help="Column names of covariates (separate by comma)"),
    make_option("--sv_output_path", help="Output file for single variant analysis (%CHR% and %PHENO% substituted)"),
    make_option("--group_output_path", help="Output file for group-based tests (%CHR% and %PHENO% substituted)"),
    make_option("--skat_o_method", help="Method for 'skatOMeta' function' ('integration'/'saddlepoint'). If lowest p-value in T1 or SKAT test < 1E-9, change this to 'saddlepoint'", default="integration")
  )
  
  parse_args(OptionParser(option_list=option_list))
}

make_filename = function(path, chr, pheno) {
  path = gsub("%CHR%", chr, path)
  gsub("%PHENO%", pheno, path)
}

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
  
  list(genotype_matrix = geno_matrix, phenotype_matrix = pheno_matrix)
}

build_snpinfo = function(gene, snps) {
  snp_info = data.frame(
    gene = as.character(rep(gene, length(snps))),
    Name = as.character(snps))
  snp_info
}

calculate_null_model_residuals = function(phenotype_matrix, model_formula) {
  null_model = stats::glm(formula = model_formula, data = phenotype_matrix)
  residuals = stats::residuals(null_model, type = "response")
  print(paste("Calculated ", length(residuals), " residuals for null model.", sep=""))
  residuals
}

process_gene = function(parameters, gene, snps, phenotype, write_header = F) {
  print(paste("==== Process gene ", gene, " (", length(snps), " variants, chr", parameters$chr, ")", sep=""))
  
  genotype = load_genotype(parameters$chr, parameters$bgen_path, snps)
  geno_pheno = prepare_genotype_phenotype_matrices(genotype, phenotype)
  snp_info = build_snpinfo(gene, snps)
  #residuals = calculate_null_model_residuals(geno_pheno$phenotype_matrix, parameters$model_formula)

  family = "gaussian"
  if (parameters$phenotype_type == "binary") {
    family = "binomial"
  }
  
  scores = prepScores2(Z = geno_pheno$genotype_matrix,
                       formula = parameters$model_formula,
                       family = family,
                       SNPInfo = snp_info,
                       data = geno_pheno$phenotype_matrix)
  results = perform_tests(scores, snp_info, parameters$skat_o_method)
  
  write_sv_result(results$single_variant, parameters$sv_output_file, write_header)
  write_group_result(results, parameters$group_output_file, write_header)
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

prepare_phenotype = function(phenotype_file, phenotype_col, individual_col, covariate_cols, phenotype_type) {
  print("======= PREPARING PHENOTYPES")
  
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

  print(paste("Phenotype type: ", phenotype_type, sep=""))
  if (phenotype_type == "binary") {
    pheno_levels = levels(as.factor(phenotypes2[,phenotype_col]))
    print(paste("Levels for binary phenotype: ", paste(pheno_levels, collapse=", ", sep=""), sep=""))
    if (length(pheno_levels) != 2) {
      stop("Require exactly two levels for binary phenotype.")
    }
    print(table(phenotypes2[,phenotype_col]))
  }

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
  
  print("SKAT test")
  skat = skatMeta(scores, SNPInfo = snp_info)
  print(format(head(skat), digits = 2))
  
  print("Burden test")
  burden = burdenMeta(scores, SNPInfo = snp_info)
  print(format(head(burden), digits=2))
  
  print("SKAT-O test")
  skat_o = skatOMeta(scores, SNPInfo=snp_info, burden.wts = function(maf) { dbeta(maf, 1, 25) }, method = skat_o_method)
  print(format(head(skat_o), digits=2))
  
  list(single_variant = single_variant, burden = burden, skat = skat, skat_o = skat_o)
}

check_column = function(headers, col_name) {
  if (!(col_name %in% colnames(headers))) {
    stop(paste("Phenotype file does not contain column:", col_name))
  }
}

check_and_prepare_parameters = function(parameters) {
  if (length(parameters$chr) == 0) {
    stop("CHR parameter missing - maybe try option '--help'")
  }
  
  parameters$sv_output_file = make_filename(parameters$sv_output_path, parameters$chr, parameters$phenotype_col)
  #if (file.access(parameters$sv_output_file, mode=2) != 0) {
  #  stop(paste("Single-variant output file not writeable (check path): ", parameters$sv_output_file, sep=""))
  #}

  parameters$group_output_file = make_filename(parameters$group_output_path, parameters$chr, parameters$phenotype_col)
  #if (file.access(parameters$group_output_file, mode=2) != 0) {
  #  stop(paste("Group output file not writeable (check path): ", parameters$group_output_file, sep=""))
  #}
  
  bgen_file = make_filename(parameters$bgen_path, parameters$chr, parameters$phenotype_col)
  if (!file.exists(bgen_file)) {
    stop(paste("BGEN file not found:", bgen_file))
  }
  
  bgen_bgi_file = paste(bgen_file, ".bgi", sep="")
  if (!file.exists(bgen_file)) {
    stop(paste("BGEN index file (bgi) not found:", bgen_bgi_file))
  }

  if (!file.exists(parameters$group_file)) {
    stop(paste("Group file not found:", parameters$group_file))
  }
  
  if (!file.exists(parameters$phenotype_file)) {
    stop(paste("Phenotype file not found:", parameters$phenotype_file))
  }

  headers = read.table(parameters$phenotype_file, h=T, sep="\t", nrows = 5)
  if (nrow(headers) != 5) {
    stop(paste("Phenotype file invalid (expected to get 5 rows):", parameters$phenotype_file))
  }
  
  check_column(headers, parameters$phenotype_col)

  if (parameters$phenotype_type != "quantitative" && parameters$phenotype_type != "binary") {
    stop("Phenotype type must be 'binary' or 'quantitative'.")
  }
  
  covariates_arr = trimws(unlist(strsplit(parameters$covariate_cols, ",", fixed=T)[[1]]))
  for (adj_col in covariates_arr) {
    check_column(headers, adj_col)
  }
  
  if (!(parameters$skat_o_method %in% c("integration", "saddlepoint"))) {
    stop(paste("Invalid skat_o_method (must be 'integration' or 'saddlepoint'):", parameters$skat_o_method))
  }

  parameters$model_formula = prepare_formula(parameters$phenotype_col, parameters$covariate_cols)

  print("Parameter dump:")
  print(parameters)  
  parameters
}

clean_previous_output = function(parameters) {
  if (file.exists(parameters$sv_output_file)) {
    print("Clean previous SV output file")
    file.remove(parameters$sv_output_file)
  }

  if (file.exists(parameters$group_output_file)) {
    print("Clean previous group output file")
    file.remove(parameters$group_output_file)
  }
}

process_group_file = function(parameters, phenotype) {
  print("======= PROCESS GROUP FILE")
  
  con = file(parameters$group_file, "r")
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
    if (parameters$chr == line_chr) {
      process_gene(parameters, gene, snps, phenotype, write_header)
      write_header = F
    }
  }
  close(con)
  
  print("Finished")
}

perform_analysis = function() {
  parameters = parse_options()
  parameters = check_and_prepare_parameters(parameters)
  clean_previous_output(parameters)
  phenotype = prepare_phenotype(parameters$phenotype_file, parameters$phenotype_col, "individual_id", parameters$covariate_cols, parameters$phenotype_type)
  process_group_file(parameters, phenotype)
}

perform_analysis()
