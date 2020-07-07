CREATE TABLE exons (
  ensembl_gene_id varchar(30),
  ensembl_chrom_start int,
  ensembl_chrom_end int,
  ensembl_exon_id varchar(30),
  rank int
);

CREATE INDEX exon_gene_id on exons (ensembl_gene_id);


