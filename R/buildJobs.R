#!/usr/local/bin/Rscript

library(optparse)

parse_options = function() {
	option_list = list(
		make_option("--parameters", help="Parameters input file path"),
		make_option("--jobs", help="Jobs output file path")
	)
	args = parse_args(OptionParser(option_list=option_list))
	return(args)
}


read_parameters = function(parameters_fn) {
	if (parameters_fn == "" || length(parameters_fn) == 0) {
		stop("Need parameters file name.")
	}
	if (!file.exists(parameters_fn)) {
		stop(paste("Parameter file not found: ", parameters_fn, sep=""))
	}
	parameters = read.table(parameters_fn, sep="=", as.is=T, h=F)
	colnames(parameters) = c("param_key", "param_value")
	print("Parameters:")
	for (i in 1:nrow(parameters)) {
		print(paste(parameters[i, "param_key"], ": ", 
			parameters[i, "param_value"], sep=""))
	}
	parameters
}

get_param = function(parameters, key, optional=F) {
	line = which(parameters$param_key == key)
	if (length(line) == 0 && !optional) {
		stop(paste("Parameter missing: ", key, sep=""))
	}
	parameters[line, "param_value"]
}

read_phenotype_file = function(parameters) {
	fn = get_param(parameters, "PHENOTYPE_FILE")
	print(paste("Reading phenotype file: ", fn, sep=""))
	if (!file.exists(fn)) {
		stop(paste("Phenotype file not found: ", fn, sep=""))
	}

	data = read.table(fn, h=T)
	if (nrow(data) == 0) {
		stop("No phenotype data.")
	}

	print(paste("Got phenotype rows: ", nrow(data), sep=""))
	data
}

parse_phenotypes = function(parameters) {
        phenotype_columns_str = get_param(parameters, "PHENOTYPES")
        strsplit(phenotype_columns_str, ",", fixed=T)[[1]]
}

check_phenotype_vars = function(parameters, data) {
	phenotype_columns = parse_phenotypes(parameters)
	if (length(phenotype_columns) == 0) {
		stop("No phenotypes given - use parameter PHENOTYPES.")
	}

	available_columns = colnames(data)

	for (phenotype in phenotype_columns) {
		print(paste("Check phenotype column: ", phenotype, sep=""))
		if (!(phenotype %in% available_columns)) {
			stop(paste("Phenotype not found: ", phenotype, sep=""))
		}
	}

	phenotype_columns
}

check_global_covariates = function(parameters, data) {
	covariables_str = get_param(parameters, "COVARIATES")
        covariables_columns = strsplit(covariables_str, ",", fixed=T)[[1]]
	if (length(covariables_columns) == 0) {
		print("No global covariates.")
		return()
	}

        available_columns = colnames(data)

        for (covariable in covariables_columns) {
                print(paste("Check covariate column: ", covariable, sep=""))
                if (!(covariable %in% available_columns)) {
                        stop(paste("Covariable not found: ", phenotype, sep=""))
                }
        }
}

check_specific_covariates_for_phenotype = function(parameters, data, phenotype) {
        covariables_str = get_param(parameters, paste("COVARIATE_", phenotype, sep=""), optional=T)
	if (length(covariables_str) == 0) {
                print(paste("No specific covariates for phenotype: ", phenotype, sep=""))
                return()
        }

        covariables_columns = strsplit(covariables_str, ",", fixed=T)[[1]]

        available_columns = colnames(data)

        for (covariable in covariables_columns) {
                print(paste("Check covariate column: ", covariable, 
			"; phenotype: ", phenotype, sep=""))
                if (!(covariable %in% available_columns)) {
                        stop(paste("Covariable not found: ", phenotype,
			"; phenotype: ", phenotype, sep=""))
                }
        }
}

check_specific_covariates = function(parameters, data, phenotype_columns) {
	for (phenotype in phenotype_columns) {
		print(paste("Checking additional covariates for: ", phenotype, sep=""))
		check_specific_covariates_for_phenotype(parameters, data, phenotype)
	}
}

check_covariates = function(parameters, data, phenotype_columns) {
	check_global_covariates(parameters, data)
	check_specific_covariates(parameters, data, phenotype_columns)
}

check_individual_id = function(parameters, data) {
	id_col = get_param(parameters, "INDIVIDUAL_ID_COLUMN")
        available_columns = colnames(data)
        if (!(id_col %in% available_columns)) {
                stop(paste("Individual ID column not found: ", id_col, sep=""))
        } else {
                print(paste("Individual ID column found: ", id_col, sep=""))
	}
}	

check_phenotypes_and_covariates = function(parameters) {
	data = read_phenotype_file(parameters)
	phenotype_columns = check_phenotype_vars(parameters, data)
	check_covariates(parameters, data, phenotype_columns)
	check_individual_id(parameters, data)
}	

parse_groups = function(parameters) {
	groups_str = get_param(parameters, "GROUPS")
        groups_arr = strsplit(groups_str, ",", fixed=T)[[1]]
	if (length(groups_arr) == 0) {
		stop(paste("Need at least one group: ", groups_str, sep=""))
	}
	groups_arr
}

check_group_file = function(parameters, group) {
	print(paste("Check group: ", group, sep=""))
	group_fn = get_param(parameters, paste("GROUP_PATH_", group, sep=""))
	print(paste("Group file: ", group_fn))
	if (!file.exists(group_fn)) {
		stop(paste("Group file not found: ", group_fn, sep=""))
	}
}

check_groups = function(parameters) {
	groups = parse_groups(parameters)
	for (group in groups) {
		check_group_file(parameters, group)
	}
}

check_path = function(parameters, param) {
	fn = get_param(parameters, param)
	if (!file.exists(fn)) {
		stop(paste("Required file missing: ", param, " / ", fn, sep=""))
	}
}

check_paths = function(parameters) {
	check_path(parameters, "GROUP_TEST_SCRIPT_PATH")
	check_path(parameters, "COLLECT_SCRIPT_PATH")
	check_path(parameters, "MART_MAPPING_FILE")
}

check_parameters = function(parameters) {
	check_groups(parameters)
	check_phenotypes_and_covariates(parameters)
	check_paths(parameters)
}

read_and_check_parameters = function(args) {
	parameters = read_parameters(args$parameters)
	check_parameters(parameters)
	return(parameters)
}

make_covar_cols = function(parameters, phenotype) {
	global_covars = get_param(parameters, "COVARIATES")
	local_covars = get_param(parameters, paste("COVARIATE_", phenotype, sep=""), optional=T)
        if (length(local_covars) == 0) {
                return(global_covars)
        } else {
		return(paste(global_covars, local_covars, sep=","))
	}
}

build_group_test_jobs = function(parameters, group, phenotype, jobs_fn) {
	script_path = get_param(parameters, "GROUP_TEST_SCRIPT_PATH")
	bgen_path = get_param(parameters, "BGEN_PATH")
	pheno_path = get_param(parameters, "PHENOTYPE_FILE")
	output_dir = get_param(parameters, "OUTPUT_DIRECTORY")
	log_dir = get_param(parameters, "LOG_DIRECTORY")

	dir.create(output_dir, showWarnings = F)
	dir.create(log_dir, showWarnings = F)

	job_name = paste(group, phenotype, sep="_")

	covar_cols = make_covar_cols(parameters, phenotype)

	group_fn = get_param(parameters, paste("GROUP_PATH_", group, sep=""))

	for (chr in 1:22) {
		log_fn = paste(log_dir, "/", job_name, "-chr", chr, "-%j.log", sep="")
		command = paste("sbatch \\\n",
			" --job-name=", job_name, " \\\n",
			" --output=\"", log_fn, "\" \\\n",
			" --error=\"", log_fn, "\" \\\n",
			" ", script_path, " \\\n",
			" --chr=", chr, " \\\n",
			" --bgen_path=\"", bgen_path, "\" \\\n",
			" --group_file=\"", group_fn, "\" \\\n",
			" --phenotype_file=\"", pheno_path, "\" \\\n",
			" --phenotype_col=", phenotype, " \\\n",
			" --covariate_cols=\"", covar_cols, "\" \\\n",
			" --sv_output_path=\"", output_dir, "/sv-", group, "-%PHENO%-chr%CHR%.txt\" \\\n",
			" --group_output_path=\"", output_dir, "/group-", group, "-%PHENO%-chr%CHR%.txt\"\n",
			"\n\n",
			sep="")
		cat(command, file = jobs_fn, append = T)
	}
}

build_collect_and_plot_job = function(parameters, group, phenotype, jobs_fn) {
	script_path = get_param(parameters, "COLLECT_SCRIPT_PATH")
	mart_mapping_file = get_param(parameters, "MART_MAPPING_FILE")
	filter_formula = get_param(parameters, "FILTER_FORMULA")

        output_dir = get_param(parameters, "OUTPUT_DIRECTORY")
        log_dir = get_param(parameters, "LOG_DIRECTORY")
        job_name = paste(group, phenotype, sep="_")
        log_fn = paste(log_dir, "/", job_name, "-collect-%j.log", sep="")

	qq_fn = paste(output_dir, "/", group, "-",
		phenotype, "-qq.pdf", sep="")

	xlsx_fn = paste(output_dir, "/", group, "-",
		phenotype, "-top.xlsx", sep="")

        group_result_fn = paste(output_dir, "/group-", group, "-", 
		phenotype, "-chr%CHR%.txt", sep="")

        command = paste("sbatch \\\n",
		" --dependency=singleton \\\n",
                " --job-name=", job_name, " \\\n",
                " --output=\"", log_fn, "\" \\\n",
                " --error=\"", log_fn, "\" \\\n",
		" ", script_path, " \\\n",
		" --group_result_files=\"", group_result_fn, "\" \\\n",
		" --mart_mapping_file=", mart_mapping_file, " \\\n",
		" --qq_plot_output_file=\"", qq_fn, "\" \\\n",
		" --top_xlsx_output_file=\"", xlsx_fn, "\" \\\n",
		" --filter_formula=\"", filter_formula, "\"\n",
		"\n\n",
		sep="")

	cat(command, file = jobs_fn, append = T)
}

build_jobs_for_phenotype = function(parameters, groups, phenotype, jobs_fn) {
	for (group in groups) {
		cat(paste("#### ", phenotype, "\n", sep=""), file=jobs_fn, append=T)
		build_group_test_jobs(parameters, group, phenotype, jobs_fn)
		build_collect_and_plot_job(parameters, group, phenotype, jobs_fn)
	}
}

write_jobs_header = function(jobs_fn) {
	cat("#!/bin/bash\n\n", file=jobs_fn)
}

build_jobs = function(parameters, jobs_fn) {
	write_jobs_header(jobs_fn)
	phenotypes = parse_phenotypes(parameters)
	groups = parse_groups(parameters)
	for (phenotype in phenotypes) {
		build_jobs_for_phenotype(parameters, groups, phenotype, jobs_fn)
	}
	print(paste("Built jobs file: ", jobs_fn, sep=""))
}

read_parameters_and_build_jobs = function() {
	args = parse_options()
	parameters = read_and_check_parameters(args)
	build_jobs(parameters, args$jobs)
}

read_parameters_and_build_jobs()
