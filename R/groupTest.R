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
    make_option("--kinship", help="Kinship table path (requires columns 'ID1', 'ID2' and 'Kinship'); default '' = don't use matrix", default="")
  )
  
  parse_args(OptionParser(option_list=option_list))
}

make_filename = function(path, chr, pheno) {
  path = gsub("%CHR%", chr, path)
  gsub("%PHENO%", pheno, path)
}

load_genotype = function(chr, bgen_path, snps, min_maf, max_maf) {
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
      colnames(chunk_geno)[2] = snps[1]
    } else {
      chunk_geno = data.frame(individual_id = data$samples)
    }
    
    if (nrow(my_geno) == 0) {
      my_geno = chunk_geno
    } else {
      if (nrow(chunk_geno) != length(data$samples)) {
        stop("Missing samples in chunk.")
      }
      if (ncol(chunk_geno) > 1) {
        my_geno = cbind(my_geno, chunk_geno[, 2:ncol(chunk_geno)])
      }
    }
  }
  
  print(paste("Genotype data frame (all chunks): ", dim(my_geno)[1], "x", dim(my_geno)[2], sep=""))

  if (ncol(my_geno) > 1) {
    if (ncol(my_geno) > 2) {
      mafs = rowMeans(sum_dosages, na.rm=T)/2
    } else {
      mafs = c(mean(sum_dosages)/2)
    }
    fail = which(mafs < min_maf | mafs > max_maf)
    print(paste(length(fail), " SNPs fail MAF range filter", sep=""))
    if (length(fail) > 0) {
      if (length(fail) + 1 == ncol(my_geno)) {
        my_geno = data.frame(individual_id = data$samples)
      } else {
        my_geno = my_geno[,-(fail+1)] # col 1 is IID
      }
      print(paste("Genotype data frame after MAF filter: ", dim(my_geno)[1], "x", dim(my_geno)[2], sep=""))
    }
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

build_kinship_matrix = function(kinship_list, genotype_matrix) {
  dimid = colnames(genotype_matrix)
  id1 = match(kinship_list$ID1, dimid)
  id2 = match(kinship_list$ID2, dimid)
  m = sparseMatrix(i=id1, j=id2, x=x, #symmetric=TRUE, 
               dims = c(length(dimid), length(dimid)),
               dimnames=list(as.character(dimid), as.character(dimid)))
  return(m)
}

process_gene = function(parameters, gene, snps, phenotype, kinship, write_header = F) {
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
                           as.numeric(parameters$max_maf))

  if (ncol(genotype) <= 2) {
    print(paste("Need more than one variant, but got:", ncol(genotype)-1))
    return(F)
  }

  geno_pheno = prepare_genotype_phenotype_matrices(genotype, phenotype)
  snp_info = build_snpinfo(gene, snps)
  #residuals = calculate_null_model_residuals(geno_pheno$phenotype_matrix, parameters$model_formula)

  family = "gaussian"
  if (parameters$phenotype_type == "binary") {
    family = "binomial"
  }
  
  if (nrow(kinship) > 0) {
    kinship_matrix = build_kinship_matrix(kinship, geno_pheno$genotype_matrix)
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
  results = perform_tests(scores, snp_info, parameters$skat_o_method,
			  as.numeric(parameters$min_maf), as.numeric(parameters$max_maf))
  
  write_sv_result(results$single_variant, parameters$sv_output_file, write_header)
  write_group_result(results, nrow(geno_pheno$phenotype_matrix), parameters$group_output_file, write_header)

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

perform_tests = function(scores, snp_info, skat_o_method, min_maf, max_maf) {
  single_variant = singlesnpMeta(scores, SNPInfo = snp_info)
  
  print("SKAT test")
  skat = skatMeta(scores, SNPInfo = snp_info)
  print(format(head(skat), digits = 2))
  
  print("Burden test")
  burden = burdenMeta(scores, SNPInfo = snp_info, mafRange = c(min_maf, max_maf))
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

process_group_file = function(parameters, phenotype, kinship) {
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
      got_output = process_gene(parameters, gene, snps, phenotype, kinship, write_header)
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
    return(kinship)
  } else {
    print("Not using a kinship matrix.")
    return(data.frame())
  }
}

perform_analysis = function() {
  parameters = parse_options()
  parameters = check_and_prepare_parameters(parameters)
  clean_previous_output(parameters)
  phenotype = prepare_phenotype(parameters$phenotype_file, parameters$phenotype_col, parameters$individual_col, parameters$covariate_cols, parameters$phenotype_type, parameters$exclude_individuals)
  kinship = prepare_kinship_file(parameters$kinship)
  process_group_file(parameters, phenotype, kinship)
}

perform_analysis()
