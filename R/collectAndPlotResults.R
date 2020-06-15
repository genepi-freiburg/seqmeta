#!/usr/local/R/R-3.6.3/bin/Rscript

print("Loading R packages")
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(qqman))
suppressPackageStartupMessages(library(xlsx))

parse_options = function() {
  option_list = list(
    make_option("--group_result_files", help="Files containing group-based test results (needs to have %CHR% if files are split by chromosome)"),
    make_option("--qq_plot_output_file", help="Path/filename of PDF output file for QQ plot", default="qq_plots.pdf"),
    make_option("--top_xlsx_output_file", help="Path/filename of XLSX output file for top results", default="top_results_burden_p0.05.xlsx"),
    make_option("--mart_mapping_file", help="File with ENSG mappings", default="mart_export.txt")
  )

  parse_args(OptionParser(option_list=option_list))
}


check_options = function(options) {
	if (!file.exists(options$mart_mapping_file)) {
		stop(paste("MART mapping file does not exist:", options$mart_mapping_file))
	}

	if (length(options$group_result_files) == 0) {
		stop("Need 'group_result_files' parameter.")
	}

	for (chr in 1:22) {
		fn = gsub("%CHR%", chr, options$group_result_files)
		if (!file.exists(fn)) {
			stop(paste("Group results file not found:", fn))
		}
	}
}

read_data = function(group_output_files) {
	data = data.frame()

	for (chr in 1:22) {
		print(paste("Collect chromosome", chr))
		fn = gsub("%CHR%", chr, group_output_files)
		d = read.table(fn, h=T)
		data = rbind(d, data)
	}

	print(summary(data))

	data
}

do_qq_plot = function(p, title) {
	print(paste("Plot QQ plot:", title, "for", length(p), "p-values"))
	p2 = p[!is.na(p)]
	print(paste("Non-NA:", length(p2)))
	print(summary(p2))

	if (length(which(p2 > 1)) > 0) {
		print("Funny p-values > 1:")
		print(p2[which(p2 > 1)])
		p2[which(p2 > 1)] = 1
	}
	
	chisq = qchisq(1 - p2, 1)
        lambda = median(chisq) / qchisq(0.5, 1)
	print(paste("lambda:", lambda))

	qq(p, main = title, pch = 18, col = "blue4", cex = 1.5, las = 1,
	   sub = paste(length(p2), "p-values; lambda =", round(lambda, 2)))
}

do_qq_plots = function(qq_plot_out_file, data) {
	pdf(qq_plot_out_file)

	do_qq_plot(data$burden_p, title = "QQ plot, Burden test")
	do_qq_plot(data$skat_p, title = "QQ plot, SKAT test")
	do_qq_plot(data$skat_o_p, title = "QQ plot, SKAT-O test")

	dev.off()
}

do_top_results = function(top_out_file, mart_mapping_file, data) {
	top_results = subset(data, data$burden_p < 0.05)

	mart = read.table(mart_mapping_file, h = T, sep="\t")
	colnames(mart) = c("gene", "Gene_Symbol", "Gene_Start_b38", "Gene_Chr")

	top_results_annotated = merge(top_results, mart, all.x=T)
	top_results_annotated = top_results_annotated[order(top_results_annotated$burden_p), ]

	if (length(which(is.na(top_results_annotated$GeneSymbol))) > 0) {
		print("There are gene IDs without symboles.")
	}

	print("Writing XLSX")
	write.xlsx(top_results_annotated, top_out_file, row.names=F, col.names=T)
}

plot_and_annotate = function() {
	options = parse_options()
	check_options(options)

	group_output_files = options$group_result_files
	qq_plot_out_file = options$qq_plot_output_file
	top_out_file = options$top_xlsx_output_file
	mart_mapping_file = options$mart_mapping_file

	data = read_data(group_output_files)
	do_qq_plots(qq_plot_out_file, data)
	do_top_results(top_out_file, mart_mapping_file, data)
}

plot_and_annotate()
