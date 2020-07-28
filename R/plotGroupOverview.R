# script to plot histograms for group files: variant count per gene

dirs = c("01_ultraRareDamaging", "02_synonymous", "03_rareDamaging")
library(ggplot2)

pdf("variants_per_gene.pdf")
for (dir in dirs) {
  file = paste(dir, "variants_and_genes.txt", sep="/")
  print(paste("Read: ", file, sep=""))
  data = read.table(file, h=F, sep="\t")
  colnames(data) = c("VARIANT", "GENE")
  print(paste("Got rows: ", nrow(data), sep=""))
  data$GENE = as.factor(data$GENE)
  print(summary(data))
  print(head(data))
  data_table = data.frame(table(data$GENE))
  print(head(data_table))
  plot = ggplot(data_table, aes(x = Freq)) + geom_histogram(color="darkblue", fill="lightblue") + scale_y_log10() + ggtitle(dir) + xlab("Variants per Gene") + ylab("Frequency")
  print(plot)
  plot
}
dev.off()

