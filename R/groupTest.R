#!/usr/local/R/R-3.6.3/bin/Rscript

# can be set to F to remove dependency to pryr package
SHOW_MEMORY_USAGE = T

# read in BGEN files in chunks of x SNPs (rbgen error)
BGEN_CHUNK_SIZE = 1000

print("Loading packages")
suppressPackageStartupMessages(library(rbgen))
suppressPackageStartupMessages(library(seqMeta))
suppressPackageStartupMessages(library(optparse))
if (SHOW_MEMORY_USAGE) {
  suppressPackageStartupMessages(library("pryr"))
}
suppressPackageStartupMessages(library(Matrix))

parse_options = function() {
  option_list = list( 
    make_option("--chr", help="Chromosome number (for patterns)"),
    make_option("--bgen_path", help="BGEN/BGI file name, %CHR% will be substituted by 'chr'"),
    make_option("--sample_path", help="Sample file name, %CHR% will be substituted by 'chr'", default=""),
    make_option("--group_file", help="Group file name"),
    make_option("--phenotype_file", help="Phenotype file name"),
    make_option("--exclude_individuals", help="File name with individuals to exclude, one in a row. Default empty = no exclusions", default=""),
    make_option("--phenotype_col", help="Column name of phenotype"),
    make_option("--individual_col", help="Column name of individual ID, default: individual_id", default="individual_id"),
    make_option("--phenotype_type", help="Phenotype type ('binary'/'quantitative'), default: quantitative", default="quantitative"),
    make_option("--covariate_cols", help="Column names of covariates (separate by comma)"),
    make_option("--sv_output_path", help="Output file for single variant analysis (%CHR% and %PHENO% substituted)"),
    make_option("--group_output_path", help="Output file for group-based tests (%CHR% and %PHENO% substituted)"),
    make_option("--skat_o_method", help="Method for 'skatOMeta' function' ('integration'/'saddlepoint'). If lowest p-value in T1 or SKAT test < 1E-9, change this to 'saddlepoint'", default="integration"),
    make_option("--min_maf", help="Lower minor allele frequency (MAF) treshhold, e.g. 0.01. Default 0, range [0,1].", default="0"),
    make_option("--max_maf", help="Upper minor allele frequency (MAF) treshhold, e.g. 0.01. Default 1, range [0,1].", default="1"),
    make_option("--kinship", help="Kinship table path (requires columns 'ID1', 'ID2' and 'Kinship'); default '' = don't use matrix", default=""),
    make_option("--stepwise_grouptest", help="Step-wise group test ('add-one-in'); give gene name to calculate", default=""),
    make_option("--step_p_limit", help="P value threshold for single variants to be included in step-wise group test; default: 0.1", default="0.1"),
    make_option("--leave_one_out", help="Use add-one-in mode (0) or leave-one-out mode (1); default: 0", default="0"),
    make_option("--skip_skat", help="Skip the SKAT test (1), default: 0", default="0"),
    make_option("--skip_skato", help="Skip the SKAT-O test (1), default: 0", default="0"),
    make_option("--export_matrices", help="Exports genotype and phenotype matrix for given gene(s).", default="")
  )
  
  parse_args(OptionParser(option_list=option_list))
}

make_filename = function(path, chr, pheno) {
  path = gsub("%CHR%", chr, path)
  gsub("%PHENO%", pheno, path)
}

load_genotype = function(chr, bgen_path, snps, min_maf, max_maf, sample) {
  bgen_file = gsub("%CHR%", chr, bgen_path)
  print(paste("Loading BGEN file: ", bgen_file, sep=""))
  
  chunk_count = ceiling(length(snps) / BGEN_CHUNK_SIZE)
  print(paste("Load ", length(snps), " SNPs in ", chunk_count, " chunks.", sep=""))

  my_geno = data.frame()
  
  for (chunk_index in 1:chunk_count) {
    startOffset = BGEN_CHUNK_SIZE * (chunk_index - 1) + 1
    stopOffset = min(startOffset + BGEN_CHUNK_SIZE - 1, length(snps))
    
    print(paste("Chunk ", chunk_index, ": SNPs ", startOffset, " to ", stopOffset, sep=""))
    data = bgen.load(bgen_file, rsids = snps[startOffset:stopOffset])
    print(paste("Got data for ", length(data$samples), " samples and ", length(data$variants$rsid), " variants", sep=""))
    
    sum_dosages = 0 * data$data[,,3] + data$data[,,2] + 2 * data$data[,,1]
    if (length(data$variants$rsid) > 1) {
      sd2 = t(sum_dosages)
      chunk_geno = cbind(data$samples, sd2)
      colnames(chunk_geno)[1] = "individual_id"
    } else if (length(data$variants$rsid) == 1) {
      chunk_geno = data.frame(individual_id = data$samples, snp = sum_dosages)
      colnames(chunk_geno)[2] = snps[startOffset]
    } else {
      chunk_geno = data.frame(individual_id = data$samples)
    }

    if (ncol(chunk_geno) > 2) {
      mafs = rowMeans(sum_dosages, na.rm=T)/2
    } else {
      mafs = c(mean(sum_dosages)/2)
    }

    fail = which(mafs < min_maf | mafs > max_maf)
    print(paste(length(fail), " SNPs fail MAF range filter", sep=""))
    if (length(fail) > 0) {
      if (length(fail) + 1 == ncol(chunk_geno)) {
        chunk_geno = data.frame(individual_id = data$samples)
      } else {
        chunk_geno = chunk_geno[,-(fail+1)] # col 1 is IID
      }
      print(paste("Chunk genotype data frame after MAF filter: ", dim(chunk_geno)[1], "x", dim(chunk_geno)[2], sep=""))
    }
    
    if (nrow(my_geno) == 0) {
      my_geno = chunk_geno
    } else {
      if (nrow(chunk_geno) != length(data$samples)) {
        stop("Missing samples in chunk.")
      }
      if (ncol(chunk_geno) > 2) {
        my_geno = cbind(my_geno, chunk_geno[, 2:ncol(chunk_geno)])
      } else if (ncol(chunk_geno) == 2) {
	chunk_geno_frame = data.frame(V1 = chunk_geno[,2])
        colnames(chunk_geno_frame) = snps[startOffset]
        print(paste("assigned singular SNP name:", snps[startOffset]))
        my_geno = cbind(my_geno, chunk_geno_frame)
      }
    }
  }
  
  print(paste("Genotype data frame (all chunks): ", dim(my_geno)[1], "x", dim(my_geno)[2], sep=""))

  if (nrow(sample) > 0) {
    print(paste("Apply sample file with", nrow(sample), "rows"))
    if (nrow(sample) != nrow(my_geno)) {
      stop("Sample file size / genotype file size mismatch")
    }
    my_geno[,1] = sample$ID2
  }
  my_geno
}

prepare_genotype_phenotype_matrices = function(genotype, phenotype) {
  data = merge(phenotype, genotype, by="individual_id")
  print(paste("Merged genotype and phenotype data: ", dim(data)[1], "x", dim(data)[2], sep=""))
  if (dim(data)[1] == 0) {
    print("Merge problems - individual ID mismatch?")
    print(head(phenotype))
    print(head(genotype))
    stop()
  }
  
  if (ncol(phenotype) + 1 < ncol(data)) {
    geno_matrix = data[, (ncol(phenotype)+1):ncol(data)]
  } else {
    geno_matrix = as.matrix(data[, (ncol(phenotype)+1):ncol(data)])
    colnames(geno_matrix)[1] = colnames(data)[ncol(phenotype)+1]
  }
  
  rownames(geno_matrix) = data$individual_id
  
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

build_kinship_matrix = function(kinship_list, genotype_ids) {
  print(paste("Got ", length(genotype_ids), " individuals for kinship matrix.", sep=""))
  sub_kinship = kinship_list[kinship_list$ID1 %in% genotype_ids & kinship_list$ID2 %in% genotype_ids,]
  print(paste("Got ", nrow(sub_kinship), " Kinship coefficients for individual pairs.", sep=""))
  if (nrow(sub_kinship) == 0) {
    print(head(genotype_ids))
    print(head(kinship_list))
  }
  id1 = match(sub_kinship$ID1, genotype_ids)
  id2 = match(sub_kinship$ID2, genotype_ids)
  m = sparseMatrix(i=id1, j=id2, x=sub_kinship$Kinship,
               dims = c(length(genotype_ids), length(genotype_ids)),
               dimnames=list(as.character(genotype_ids), as.character(genotype_ids)))
  print(paste("Dimension of resulting sparse Kinship matrix:"))
  print(dim(m))
  return(m)
}

step_wise_group_test = function(parameters, results, geno_pheno, family, snp_info) {
  sv = results$single_variant
  sv = sv[order(sv$p),]
  sv = sv[!is.na(sv$p),]
  print(sv)

  steps = data.frame()
  for (nVar in 1:nrow(sv)) {
    if (is.na(sv[nVar, "p"]) || sv[nVar, "p"] > as.numeric(parameters$step_p_limit)) {
      print(paste("Break - p NA or p >", parameters$step_p_limit))
      break
    }

    steps[nVar, "nVar"] = nVar

    leave_one_out = as.numeric(parameters$leave_one_out) == 1

    if (!leave_one_out) {
      print(paste("add-one-in: use 1..", nVar, sep=""))
      vars = sv[1:nVar, "Name"]
    } else {
      print(paste("leave-one-out: use ", nVar, "..n", sep=""))
      vars = sv[nVar:nrow(sv), "Name"]
    }

    lastVar = sv[nVar,]
    print(paste("Use variants for step ", nVar, ": ", paste(vars, collapse=", "), sep=""))

    my_idx = which(colnames(geno_pheno$genotype_matrix) %in% vars)
    print(paste("Indices:", paste(my_idx, collapse=", ")))

    my_snpinfo = snp_info[snp_info$Name %in% vars,]
    print(my_snpinfo)

    if (length(my_idx) == 1) {
      my_geno = matrix(geno_pheno$genotype_matrix[, my_idx])
    } else {
      my_geno = geno_pheno$genotype_matrix[, my_idx]
    }
    print(dim(my_geno))
    rownames(my_geno) = rownames(geno_pheno$phenotype_matrix)
    colnames(my_geno) = my_snpinfo$Name

    print(colnames(my_geno))

    scores = prepScores2(Z = my_geno,
                         formula = parameters$model_formula,
                         family = family,
                         SNPInfo = snp_info,
                         data = geno_pheno$phenotype_matrix)
    step_results = perform_tests(parameters, scores, snp_info, parameters$skat_o_method,
                                 as.numeric(parameters$min_maf), as.numeric(parameters$max_maf))
    steps[nVar, "p_burden"] = step_results$burden$p 
    steps[nVar, "beta_burden"] = step_results$burden$beta
    steps[nVar, "se_burden"] = step_results$burden$se
    steps[nVar, "p_skat"] = step_results$skat$p
    steps[nVar, "p_skat_o"] = step_results$skat_o$p
#    steps[nVar, "p_sv"] = step_results$skat_o$p
	steps[nVar, "p_sv"] = lastVar$p
    steps[nVar, "beta_sv"] = lastVar$beta
    steps[nVar, "se_sv"] = lastVar$se
    steps[nVar, "maf_sv"] = lastVar$maf
    print(steps)
  }

  steps_fn = gsub("group-", paste("steps-", parameters$stepwise_grouptest, sep=""), parameters$group_output_file)
  print(paste("Writing to:", steps_fn))

  write.table(steps, steps_fn, row.names=F, col.names=T, sep="\t", quote=F)
}

export_matrices = function(parameters, gene, geno_pheno) {
  geno_fn = gsub("group-", paste("geno-", gene, "-", sep=""), parameters$group_output_file)
  pheno_fn = gsub("group-", paste("pheno-", gene, "-", sep=""), parameters$group_output_file)
  print(paste("Writing genotypes to: ", geno_fn, "; phenotypes to: ", pheno_fn, sep=""))

  rownames(geno_pheno$phenotype_matrix) = rownames(geno_pheno$genotype_matrix)
  write.table(geno_pheno$genotype_matrix, geno_fn, row.names=T, col.names=T, sep="\t", quote=F)
  write.table(geno_pheno$phenotype_matrix, pheno_fn, row.names=T, col.names=T, sep="\t", quote=F)
}


process_gene = function(parameters, gene, snps, phenotype, kinship, sample, write_header = F) {
  step_wise = parameters$stepwise_grouptest
  if (nchar(step_wise) > 0) {
    if (gene != step_wise) {
      print(paste("Step-wise group test mode, skip gene: ", gene, sep=""))
      return(F)
    }
  }

  export = parameters$export_matrices
  if (nchar(export) > 0) {
    genes = unlist(strsplit(export, ","))
    if (length(which(genes == gene)) == 0) {
      print(paste("Export matrices mode, skip gene: ", gene, sep=""))
      return(F)
    }
  }

  print(paste("==== Process gene ", gene, " (", length(snps), " variants, chr", parameters$chr, ")", sep=""))
  gc()
  if (SHOW_MEMORY_USAGE) {
    print(paste("Memory:", mem_used()))
  }

  if (length(snps) == 0) {
    print("No variants - return.")
    return(F)
  }
  
  genotype = load_genotype(parameters$chr, parameters$bgen_path, snps, as.numeric(parameters$min_maf),
                           as.numeric(parameters$max_maf), sample)

  if (ncol(genotype) <= 2) {
    print(paste("Need more than one variant, but got:", ncol(genotype)-1))
    return(F)
  }

  geno_pheno = prepare_genotype_phenotype_matrices(genotype, phenotype)

  if (nchar(export) > 0) {
    export_matrices(parameters, gene, geno_pheno)
    return(F)
  }

  snp_info = build_snpinfo(gene, snps)
  #residuals = calculate_null_model_residuals(geno_pheno$phenotype_matrix, parameters$model_formula)

  family = "gaussian"
  if (parameters$phenotype_type == "binary") {
    family = "binomial"
  }
  
  if (nrow(kinship) > 0) {
    kinship_matrix = build_kinship_matrix(kinship, rownames(geno_pheno$genotype_matrix))
    scores = prepScores2(Z = geno_pheno$genotype_matrix,
                         formula = parameters$model_formula,
                         family = family,
                         kins = kinship_matrix,
                         SNPInfo = snp_info,
                         data = geno_pheno$phenotype_matrix)
  } else {
    scores = prepScores2(Z = geno_pheno$genotype_matrix,
                         formula = parameters$model_formula,
                         family = family,
                         SNPInfo = snp_info,
                         data = geno_pheno$phenotype_matrix)
  }
  results = perform_tests(parameters, scores, snp_info, parameters$skat_o_method,
			  as.numeric(parameters$min_maf), as.numeric(parameters$max_maf))
  
  write_sv_result(results$single_variant, parameters$sv_output_file, write_header)
  write_group_result(results, nrow(geno_pheno$phenotype_matrix), parameters$group_output_file, write_header)

  if (nchar(step_wise) > 0) {
    print("Perform step-wise group test")
    step_wise_group_test(parameters, results, geno_pheno, family, snp_info)
  }

  if (SHOW_MEMORY_USAGE) {
    print(paste("Memory after calculating gene:", mem_used()))
  }

  return(T)
}

write_sv_result = function(sv, sv_output_file, write_header) {
  write.table(sv, sv_output_file, append = T, row.names = F, col.names = write_header, sep = "\t", quote = F)
}

write_group_result = function(result, indiv_count, group_output_file, write_header) {
  line = data.frame(gene = result$skat$gene,
                    sample_size = indiv_count,
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

prepare_phenotype = function(phenotype_file, phenotype_col, individual_col, covariate_cols, phenotype_type, exclude_file) {
  print("======= PREPARING PHENOTYPES")
  
  phenotypes = read.table(phenotype_file, h=T, sep="\t")
  print(paste("Got phenotypes for ", nrow(phenotypes), " individuals.", sep=""))
  
  pheno_na = which(is.na(phenotypes[, phenotype_col]))
  print(paste("Excluding ", length(pheno_na), " rows because of missing values for: ", phenotype_col, sep=""))

  if (length(pheno_na) > 0) {
    phenotypes = phenotypes[-pheno_na,]
  }
  
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

  if (nchar(exclude_file) > 0) {
    print(paste("Read individual exclusions: ", exclude_file, sep=""))
    exclude = read.table(exclude_file, h=F)
    print(paste("Got file with ", nrow(exclude), " rows and ", ncol(exclude), " columns.", sep=""))
    if (ncol(exclude) != 1) {
      stop("Require single column with individual ID in exclusion file.")
    }
    print(head(exclude[,1]))
    print(head(phenotypes[,individual_col]))
    excl_idx = which(phenotypes[, individual_col] %in% exclude[,1])
    if (length(excl_idx) > 0) {
      print(paste("Excluding ", length(excl_idx), " additional rows because of presence in exclusion file.", sep=""))
      phenotypes = phenotypes[-excl_idx,]
    } else {
      print("No additional exclusions based on exclusion file.")
    }
  }
  
  print(paste("Combine", individual_col, "and", phenotype_col))
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

perform_tests = function(parameters, scores, snp_info, skat_o_method, min_maf, max_maf) {
  single_variant = singlesnpMeta(scores, SNPInfo = snp_info)
  
  print("Burden test")
  burden = burdenMeta(scores, SNPInfo = snp_info, mafRange = c(min_maf, max_maf))
  print(format(head(burden), digits=2))
  
  if (parameters$skip_skat == "0") {
    print("SKAT test")
    skat = skatMeta(scores, SNPInfo = snp_info)
    print(format(head(skat), digits = 2))
  } else {
    print("Skipping SKAT test")
    skat = data.frame(gene = burden$gene, p = NA, Qmeta = NA, cmaf = NA, nmiss=NA, nsnps = NA)
  }

  if (parameters$skip_skato == "0") {
    print("SKAT-O test")
    skat_o = try(skatOMeta(scores, SNPInfo=snp_info, burden.wts = function(maf) { dbeta(maf, 1, 25) }, method = skat_o_method))
    if (is.character(skat_o)) {
      print(paste("Detected SKAT-O error:", skat_o))
      skat_o = data.frame(p = NA, pmin = NA, rho = NA, cmaf = NA, nmiss = NA, nsnps = NA, errflag = as.character(skat_o))
    }
    print(format(head(skat_o), digits=2))
  } else {
    print("Skipping SKAT-O test")
    skat_o = data.frame(p = NA, pmin = NA, rho = NA, cmaf = NA, nmiss = NA, nsnps = NA, errflag = "skipped")
  }
  
  list(single_variant = single_variant, burden = burden, skat = skat, skat_o = skat_o)
}

check_column = function(headers, col_name) {
  if (!(col_name %in% colnames(headers))) {
    print(paste("Available columns: ", colnames(headers), collapse=",", sep=""))
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

  if (parameters$min_maf < 0 || parameters$min_maf > 1) {
    stop(paste("Invalid min_maf: must be in [0, 1]!"))
  }

  if (parameters$max_maf < 0 || parameters$max_maf > 1) {
    stop(paste("Invalid max_maf: must be in [0, 1]!"))
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

process_group_file = function(parameters, phenotype, kinship, sample) {
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
      got_output = process_gene(parameters, gene, snps, phenotype, kinship, sample, write_header)
      if (got_output) {
        write_header = F
      }
    }
  }
  close(con)
  
  print("Finished")
}

prepare_kinship_file = function(kinship_path) {
  if (nchar(kinship_path) > 0) {
    kinship = read.table(kinship_path, h=T)
    print(paste("Got kinship table with", nrow(kinship), "rows."))
    print(head(kinship))
    if (!all(c("ID1", "ID2", "Kinship") %in% colnames(kinship))) {
      stop("Kinship tables needs to have columns: ID1, ID2, Kinshop")
    }
    return(kinship)
  } else {
    print("Not using a kinship matrix.")
    return(data.frame())
  }
}

load_sample_file = function(sample_file_name, chr) {
  print("======= LOADING SAMPLE FILE")
  if (nchar(sample_file_name) > 0) {
    sample = read.table(gsub("%CHR%", chr, sample_file_name), skip=2)
    print(paste("Loaded sample file with", nrow(sample), "sample IDs."))
    colnames(sample)[1:2] = c("ID1", "ID2")
    print(head(sample))
    return(sample)
  } else {
    print("Not using a sample file.")
    return(data.frame())
  }
}

perform_analysis = function() {
  parameters = parse_options()
  parameters = check_and_prepare_parameters(parameters)
  clean_previous_output(parameters)
  phenotype = prepare_phenotype(parameters$phenotype_file, parameters$phenotype_col, parameters$individual_col, parameters$covariate_cols, parameters$phenotype_type, parameters$exclude_individuals)
  sample = load_sample_file(parameters$sample_path, parameters$chr)
  kinship = prepare_kinship_file(parameters$kinship)
  process_group_file(parameters, phenotype, kinship, sample)
}

perform_analysis()
