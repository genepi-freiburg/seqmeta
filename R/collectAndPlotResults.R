library(qqman)
group_output_files = "output/group_tests_chr%CHR%.txt"
qq_plot_out_file = "qq_plots.pdf"
top_out_file = "top_results_burden_p0.05.xlsx"
mart_mapping_file = "mart_export.txt"

data = data.frame()

for (chr in 1:22) {
	print(paste("Collect chromosome", chr))
	fn = gsub("%CHR%", chr, group_output_files)
	d = read.table(fn, h=T)
	data = rbind(d, data)
}

summary(data)

pdf(qq_plot_out_file)

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

do_qq_plot(data$burden_p, title = "QQ plot, eGFR_crea, Burden test")
do_qq_plot(data$skat_p, title = "QQ plot, eGFR_crea, SKAT test")
do_qq_plot(data$skat_o_p, title = "QQ plot, eGFR_crea, SKAT-O test")

dev.off()

top_results = subset(data, data$burden_p < 0.05)

mart = read.table(mart_mapping_file, h = T, sep="\t")
colnames(mart) = c("gene", "GeneSymbol")

top_results_annotated = merge(top_results, mart, all.x=T)
top_results_annotated = top_results_annotated[order(top_results_annotated$burden_p), ]

if (length(which(is.na(top_results_annotated$GeneSymbol))) > 0) {
	print("There are gene IDs without symboles.")
}

print("Writing XLSX")
library(xlsx)
write.xlsx(top_results_annotated, top_out_file, row.names=F, col.names=T)

#write.table(top_results_annotated, top_out_file, row.names=F, col.names=T, sep="\t")

