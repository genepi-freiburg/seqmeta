# SeqMeta Pipeline

The pipeline consists of R and Bash scripts.
It depends on the following R packages: seqMeta, rbgen

## Analysis Workflow

![Workflow Figure](https://github.com/genepi-freiburg/seqmeta/blob/master/docs/workflow.PNG?raw=true)

## Group Test R Script

Most important script: groupTest.R

R scripts are executables and can be used with "--help" to show command options.


```
Usage: ./groupTest.R [options]


Options:
        --chr=CHR
                Chromosome number (for patterns)

        --bgen_path=BGEN_PATH
                BGEN/BGI file name, %CHR% will be substituted by 'chr'

        --group_file=GROUP_FILE
                Group file name

        --phenotype_file=PHENOTYPE_FILE
                Phenotype file name

        --phenotype_col=PHENOTYPE_COL
                Column name of phenotype

        --phenotype_type=PHENOTYPE_TYPE
                Phenotype type ('binary'/'quantitative'), default: quantitative

        --covariate_cols=COVARIATE_COLS
                Column names of covariates (separate by comma)

        --sv_output_path=SV_OUTPUT_PATH
                Output file for single variant analysis (%CHR% and %PHENO% substituted)

        --group_output_path=GROUP_OUTPUT_PATH
                Output file for group-based tests (%CHR% and %PHENO% substituted)

        --skat_o_method=SKAT_O_METHOD
                Method for 'skatOMeta' function' ('integration'/'saddlepoint'). If lowest p-value in T1 or SKAT test < 1E-9, change this to 'saddlepoint'

        -h, --help
                Show this help message and exit
```

## Collect and Plot script

This script combines the output of the chromosome-wise results of the groupTest.R script.
It collects all results per gene, calculates the "lambda" factor for each of the three tests, 
produces a QQ plot PDF with one page per test and writes XLSX (Excel) and TSV (tab-separated values) top result files.
These files contain all the results that satisfy a given condition (such as 'burden_p < 0.05').

Fields of the result file (can be used in the filtering condition):
* gene - Ensembl Gene ID
* skat_p, skat_Qmeta, skat_cmaf, skat_nmiss, skat_nsnps - results of the SKAT test
* burden_p, burden_beta, burden_se, burden_cmafTotal, burden_cmafUsed, burden_nsnpsTotal, burden_nsnpsUsed, burden_nmiss - results of the Burden test
* skat_o_p, skat_o_pmin, skat_o_rho, skat_o_cmaf, skat_o_nmiss, skat_o_nsnps, skat_o_errflag - results of the SKAT-O test
* Gene_Symbol - Gene symbol
* Gene_Start_b38 - Build 38 start coordinate of the gene
* Gene_Chr - chromosome of the gene

For details of the Burden, SKAT, SKAT-O fields see here: https://cran.r-project.org/web/packages/seqMeta/seqMeta.pdf

The script needs a gene mapping file mapping Ensembl Gene IDs to symbols, chromosome number and build 38 positions.
Scripts are provided to build such a mapping file using the Biomart.

```
Usage: R/collectAndPlotResults.R [options]


Options:
        --group_result_files=GROUP_RESULT_FILES
                Files containing group-based test results (needs to have %CHR% if files are split by chromosome)

        --qq_plot_output_file=QQ_PLOT_OUTPUT_FILE
                Path/filename of PDF output file for QQ plot

        --top_tsv_output_file=TOP_TSV_OUTPUT_FILE
                Path/filename of TSV output file for top results

        --top_xlsx_output_file=TOP_XLSX_OUTPUT_FILE
                Path/filename of XLSX output file for top results

        --filter_formula=FILTER_FORMULA
                Formula to filter results, default: 'burden_p < 0.05'

        --mart_mapping_file=MART_MAPPING_FILE
                File with ENSG mappings

        -h, --help
                Show this help message and exit
```

## plotGene.R script

This script produces a plot of a given gene that shows the single variants in the gene
by effect size and association p-value. It also plots the "exonic" regions.

It can be used to plot multiple genes by
* enumerating these genes on the command line, or by
* giving a top results file along with a condition on this file to select genes to plot.

The script needs to have a SQLite3 database of exon positions.
We provide a script to build such a database using the Ensembl Biomart.

![Gene Plot Figure](https://github.com/genepi-freiburg/seqmeta/blob/master/docs/genePlot.PNG?raw=true)

```
Usage: R/plotGene.R [options]


Options:
        --title=TITLE
                Plot title (Phenotype, group file)

        --genes_to_plot=GENES_TO_PLOT
                Gene(s) to plot (may use symbol or Ensembl ID); may give multiple genes separated by comma

        --top_file=TOP_FILE
                Gene(s) to plot; to be loaded from 'top results' file

        --top_file_formula=TOP_FILE_FORMULA
                Filtering formula for 'top results' file

        --sv_path=SV_PATH
                Path/filename of single variant association results (may use %CHR%)

        --mart_mapping_file=MART_MAPPING_FILE
                File with ENSG mappings

        --exon_db=EXON_DB
                SQLite3 exon DB

        --pdf_output_path=PDF_OUTPUT_PATH
                Path/filename of PDF output file. May use %SYMBOL%

        -h, --help
                Show this help message and exit
```

## buildJobs.R script

This script can be used to build a job file for Slurm or other scheduling systems.
It manages order and dependencies of the individual job steps and provides sensible parameters to the individual scripts.
The script is controlled by a "parameters" file in key-value format.

Parameter Name | Description | Example
------------ | ------------- | -------------
PHENOTYPE_FILE | phenotype file (tab-separated) | phenotype.txt
PHENOTYPES | names of the phenotypes (columns) to examine | egfr_ckdepi_creat,ckd
PHENOTYPE_TYPES | phenotype types (Q=quantitative, B=binary); one per phenotype | Q,B
INDIVIDUAL_ID_COLUMN |individual ID column name | individual_id
COVARIATES | covariates to use for all phenotypes | age_crea_serum,sex_male,U1,U2
COVARIATE_egfr_ckdepi_creat | covariates to use only for specific phenotypes | U3
GROUPS | group identifiers (shortcut) | rareDmg,syn
GROUP_PATH_rareDmg | path to group file | ...../GROUP/03_rareDamaging/group_file.txt
GROUP_PATH_syn | .../GROUP/02_synonymous/group_file.txt
FILTER_FORMULA | criteria for selecting top results (xlsx) | burden_p < 0.01
PLOT_FORMULA | criteria for selecting genes to plot | burden_nsnpsTotal > 5 & burden_p < 1e-3
BGEN_PATH | genotype path, use %CHR% | ..../ukb_wes_efe-chr%CHR%.bgen
GROUP_TEST_SCRIPT_PATH | group test script path (relative or absolute) | ../PIPELINE/R/groupTest.R
COLLECT_SCRIPT_PATH | collect and plot script path (relative or absolute) | ../PIPELINE/R/collectAndPlotResults.R
PLOT_SCRIPT_PATH | plot gene script path (relative or absolute) | ../PIPELINE/R/plotGene.R
MART_MAPPING_FILE | MART mapping file (with Ensembl Gene IDs and Gene symbols) | .../PIPELINE/biomart/mart_export_genes.txt
EXON_DB_FILE | Exon SQLite database | ..../PIPELINE/biomart/ensembl_exons.sqlite
OUTPUT_DIRECTORY | output directory, may be relative or absolute path, must not end with / | output
LOG_DIRECTORY | log directory, may be relative or absolute path, must not end with / | logs

```
Usage: R/buildJobs.R [options]


Options:
        --parameters=PARAMETERS
                Parameters input file path

        --jobs=JOBS
                Jobs output file path

        -h, --help
                Show this help message and exit
```

## BioMart helper files

### Genes Mapping File

Use the script `biomart/export_mart_genes.sh`.

This accesses: http://www.ensembl.org/biomart/martservice

This is the query:
```xml
<Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >       
 <Dataset name = "hsapiens_gene_ensembl" interface = "default" >         
  <Attribute name = "ensembl_gene_id" />          
  <Attribute name = "external_gene_name" />    
  <Attribute name = "start_position" />            
  <Attribute name = "chromosome_name" />  
 </Dataset>
</Query>'
```

### Exon position database

Use the script `biomart/export_mart_exons.sh` to download the exons file from Biomart.
Then use `biomart/init_db.sh` to process the downloaded text file and create a SQLite3 DB.

This DB uses the following schema:
```sql
CREATE TABLE exons (
  ensembl_gene_id varchar(30),
  exon_chrom_start int,
  exon_chrom_end int,
  ensembl_exon_id varchar(30),
  rank int
);

CREATE INDEX exon_gene_id on exons (ensembl_gene_id);
```

We are happy to share the resulting DB on request.
