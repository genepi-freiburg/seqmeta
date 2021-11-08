ReadMe - plotGene.R

Call: example

Rscript plotGene.R \
  --title="Disruptive missense variants" \
  --sv_path="/data/meta_analyses/11_TopMed_urate/GenePlots/Skript_MWuttke/STAAR_disruptive_missense.txt"  \
  --mart_mapping_file="/data/studies/06_UKBB/02_Projects/01_20210118_CKDGen_Kidney_Function/99_OLD/Exome_50k_Analyses/06_seqMeta/PIPELINE/biomart/mart_export_genes.txt" \
  --mysql_exon_db=wuttke \
  --mysql_exon_user=wuttke \
  --mysql_exon_password="He5oSuuquo" \
  --mysql_exon_host=biom131 \
  --pdf_output_path="output/Plot_STAAR_disruptive_missense_%SYMBOL%.pdf" \
  --genes_to_plot="ENSG00000109667,ENSG00000197891"


Output: pdf

Per gene, 2 presentations are provided. See examplarily output file.
Background: There might be situation in which the one or the other plot is more suitable, sometimes there is not much difference.

(A) including all variants and full gene presentation (unrestricted)

(B) including all variants but only relevant part of gene (restricted)
Here, the gene could be trimmed. On the respective side the gene was trimmed, “…” are added to indicate that the gene does not end there. 

In each plot is presented:

for variants:
effect estimate and SE is presented; the size of the bubble reflects the respective association p-value

for gene:
exons are presented as boxes; sometimes exons overlap and can thus not be differentiated
the arrow on top of the gene presentation indicates direction of transcription 
