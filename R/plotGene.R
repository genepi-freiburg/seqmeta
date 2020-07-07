#!/usr/local/bin/Rscript
##!/usr/local/R/R-3.6.3/bin/Rscript

print("Load R packages")
library(optparse)
#library(biomaRt)
library(DBI)

parse_options = function() {
  option_list = list(
    make_option("--title", help="Plot title (Phenotype, group file)"),
    make_option("--gene_to_plot", help="Gene to plot (may use symbol or Ensembl ID)"),
    make_option("--sv_path", help="Path/filename of single variant association results (may use %CHR%)"),
    make_option("--mart_mapping_file", help="File with ENSG mappings", default="mart_export.txt"),
    make_option("--exon_db", help="SQLite3 exon DB", default="ensembl_exons.sqlite"),
    make_option("--pdf_output_path", help="Path/filename of PDF output file. May use %SYMBOL%", default="plot_%SYMBOL%.pdf")
  )

  parse_args(OptionParser(option_list=option_list))
}

do_plot = function(opts) {

gene_to_plot = opts$gene_to_plot
if (nchar(gene_to_plot) == 0) {
  stop("Need gene")
}

sv_path = opts$sv_path
if (nchar(sv_path) == 0) {
  stop("Need path to SV association results")
}

mapping_path = opts$mart_mapping_file
if (!file.exists(mapping_path)) {
  stop("Need MART mapping file")
}

out_fn = opts$pdf_output_path

print("====== Read MART mapping table")
mapping = read.table(mapping_path, h=T, sep="\t")
colnames(mapping) = c("gene", "symbol", "gene_pos_b38", "gene_chr")
mapping$gene = as.factor(mapping$gene)
mapping$symbol = as.factor(mapping$symbol)
mapping$gene_chr = as.factor(mapping$gene_chr)
print(paste("Got rows:", nrow(mapping)))

print(paste("Lookup gene:", gene_to_plot))
mapping2 = mapping[mapping$symbol == gene_to_plot,]
if (nrow(mapping2) == 0) {
  print("Didn't find symbol, try Ensembl ID")
  mapping2 = mapping[mapping$gene == gene_to_plot,]
}
mapping = mapping2
if (nrow(mapping) != 1) {
  stop(paste("Given gene not found/not unique. Got rows:", nrow(mapping)))
}
print(mapping)

print("====== Read single-variant file")

chr = mapping$gene_chr[1]
gene = mapping$gene[1]
symbol = mapping$symbol[1]

sv_fn = gsub("%CHR%", chr, sv_path)

print(paste("Using Ensembl gene:", gene))
print(paste("Symbol:", symbol))
print(paste("Filename:", sv_fn))

print("Subset and remove missing betas")
d = read.table(sv_fn, h=T)
d$gene = as.factor(d$gene)
print(paste("All rows:", nrow(d)))
e = d[as.character(d$gene) == as.character(gene),]
print(paste("Rows for gene:", nrow(e)))
e_nonmiss = e[!is.na(e$beta),]
print(paste("Without missing betas:", nrow(e_nonmiss)))

e_nonmiss = merge(e_nonmiss, mapping)
print(paste("Merged to mapping:", nrow(e_nonmiss)))

print("Determine variant positions")
#11:47274694:G:A
for (i in 1:nrow(e_nonmiss)) {
  variant = as.character(e_nonmiss[i, "Name"])
  items = unlist(strsplit(variant, ":", fixed=T)[[1]])
  e_nonmiss[i, "chr"] = items[1]
  e_nonmiss[i, "pos"] = items[2]
  e_nonmiss[i, "all1"] = items[3]
  e_nonmiss[i, "all2"] = items[4]
}

e_nonmiss$chr = as.factor(e_nonmiss$chr)
e_nonmiss$pos = as.numeric(e_nonmiss$pos)
e_nonmiss$all1 = as.factor(e_nonmiss$all1)
e_nonmiss$all2 = as.factor(e_nonmiss$all2)

e_nonmiss$log_p = -log10(e_nonmiss$p)

print(summary(e_nonmiss))

print("====== Query SQLite3 for exons")

#ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "uswest.ensembl.org")
#exons = getBM(attributes = c('ensembl_exon_id', 'exon_chrom_start', 'exon_chrom_end', 'rank'),
#              filters = 'ensembl_gene_id', values = gene, mart = ensembl)


exon_db = dbConnect(RSQLite::SQLite(), opts$exon_db)
exons = dbGetQuery(exon_db, paste("SELECT DISTINCT exon_chrom_start, exon_chrom_end FROM exons WHERE ensembl_gene_id = \"", gene, "\"", sep=""))
dbDisconnect(exon_db)

print(paste("Got ", nrow(exons), " exon positions for ", gene, ".", sep = ""))

print("====== Plot")

out_fn = gsub("%SYMBOL%", symbol, out_fn)
print(paste("Write plot to:", out_fn))
pdf(out_fn)

x = e_nonmiss$pos/1000
y = e_nonmiss$beta

minx = round(min(c(x, exons$exon_chrom_start/1000)))
maxx = round(max(c(x, exons$exon_chrom_end/1000)))

miny = min(y - e_nonmiss$se)
maxy = max(y + e_nonmiss$se)

circle_size = sqrt(e_nonmiss$log_p/pi) # surface-proportional

symbols(x = x, y = y, 
        circles = circle_size,
        inches = 1/8,
        ann = F, 
        xlim = c(minx, maxx),
        ylim = c(miny, maxy),
        bty = "n")

arrows(x, y - e_nonmiss$se, x, y + e_nonmiss$se, length=0.05, angle=90, code=3, col="gray")

symbols(x = x, y = y, 
        circles = circle_size,
        inches = 1/8,
        ann = F, 
        bg = "steelblue2", 
        fg = "black",
        add = T,
        bty = "n")

maxMinusLog10p = max(e_nonmiss$log_p)
minMinusLog10p = min(e_nonmiss$log_p)
maxEffect = max(e_nonmiss$beta)

# e_nonmiss[1,"gene"],
title(main = paste(opts$title, e_nonmiss[1,"symbol"], sep=" / "),
      xlab = "Position [kb]", 
      ylab = "Effect Size",
      sub = paste("max(-log10(p)) = ", round(maxMinusLog10p, 1),
                  "; n = ", nrow(e_nonmiss), " SNVs", sep=""))

abline(h=0, col="darkgray")

rectHeight = 0.05 * maxEffect

rect(border = "black", col = "gray",
     xleft = exons$exon_chrom_start / 1000, xright = exons$exon_chrom_end / 1000,
     ytop = rectHeight, ybottom = -rectHeight)

legend_count = 4
legend_entries = (max(circle_size) / legend_count) * 1:legend_count
legend_labels = round((legend_entries ^ 2) * pi, 1)
  
legend(x = "topright",
       ncol = 1,
       pch = 21, 
       pt.cex = legend_entries,
       legend = legend_labels,
       pt.bg = "steelblue2",
       bg = "transparent",
       title = as.expression(bquote(paste('-log'['10']*'(p)')))
)

dev.off()

print("Finished")

}

opts = parse_options()
do_plot(opts)
