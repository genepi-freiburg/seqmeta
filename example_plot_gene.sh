mkdir -p BUN_plots
for GENE in `cat interesting_BUN_genes.txt`
do
../PIPELINE/R/plotGene.R \
        --gene_to_plot=$GENE \
        --sv_path=../20200612_BUN/output/bun_serum-sv-chr%CHR%.txt \
        --mart_mapping_file=../PIPELINE/mart_export.txt \
	--pdf_output_path=BUN_plots/BUN-plot-%SYMBOL%.pdf
done

ls -1 BUN_plots | wc -l
wc -l interesting_BUN_genes.txt

