# SeqMeta Pipeline

[![DOI](https://zenodo.org/badge/271744161.svg)](https://zenodo.org/badge/latestdoi/271744161)

The pipeline consists mainly of R scripts with a few supporting Bash scripts.
All R scripts are "executables" and can be used with `--help` to show command options.

It depends on the following R packages:
* `seqMeta` - for performing the actual group tests
* `rbgen` - for reading `bgen` files
* `optparse` - for parsing command line arguments
* `DBI` - for database access (package should come with standard R)
* `RMySQL` - for accessing MySQL databases (must be installed only if you plan to use this DB)
* `RSQLite` - for accessing SQLite3 databases (should be installed already; must be installed only if you plan to use this DB)

## Analysis Workflow

![Workflow Figure](https://github.com/genepi-freiburg/seqmeta/blob/master/docs/workflow.PNG?raw=true)

The `buildJobs.R` script can drive this process for multiple phenotypes and group files
by producing a shell script that is ready to be used with job schedulers.

## Group Test Script

The `groupTest.R` script is the most important script of the pipeline.
It reads and merges phenotypes and genotypes, taking care of missingness,
performs association tests and collects the results.
In one pass, both single-variant and group-based tests are performed.

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

        --min_maf=MIN_MAF
                Lower minor allele frequency (MAF) treshhold, e.g. 0.01. Default 0, range [0,1].

        --max_maf=MAX_MAF
                Upper minor allele frequency (MAF) treshhold, e.g. 0.01. Default 1, range [0,1].

        --kinship=KINSHIP
                Kinship table path (requires columns 'ID1', 'ID2' and 'Kinship'); default '' = don't use matrix

        --stepwise_grouptest=STEPWISE_GROUPTEST
                Step-wise group test ('add-one-in'); give gene name to calculate

        --step_p_limit=STEP_P_LIMIT
                P value threshold for single variants to be included in step-wise group test; default: 0.1

        --leave_one_out=LEAVE_ONE_OUT
                Use add-one-in mode (0) or leave-one-out mode (1); default: 0

        --export_matrices=EXPORT_MATRICES
                Exports genotype and phenotype matrix for given gene(s).

        -h, --help
                Show this help message and exit
```

## Collect and Plot results script

The `collectAndPlotResults.R` script combines the output of the chromosome-wise results of the `groupTest.R` script.
It collects all results per gene, calculates the "lambda" factor for each of the three tests, 
produces a QQ plot PDF with one page per test and writes XLSX (Excel) and TSV (tab-separated values) top result files.
These files contain all the results that satisfy a given condition (such as `burden_p < 0.05`).

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

## Plot Gene script

The `plotGene.R` script produces a plot of a given gene that shows the single variants in the gene
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
                SQLite3 exon DB, file path; required if no MySQL db given

        --mysql_exon_db=MYSQL_EXON_DB
                MySQL exon DB, database name; required if no SQLite3 DB given

        --mysql_exon_user=MYSQL_EXON_USER
                MySQL exon DB, user name; required if no SQLite3 DB given

        --mysql_exon_password=MYSQL_EXON_PASSWORD
                MySQL exon DB, password; required if no SQLite3 DB given

        --mysql_exon_host=MYSQL_EXON_HOST
                MySQL exon DB, server host; required if no SQLite3 DB given

        --pdf_output_path=PDF_OUTPUT_PATH
                Path/filename of PDF output file. May use %SYMBOL%

        -h, --help
                Show this help message and exit
```

The "sv" (single variant) file must have at least the following columns:

* gene - Ensembl Gene ID
* Name - variant name (chr:pos:ref:alt), is split to get chr/pos/ref/alt columns
* p - P-value
* beta - Effect size
* se - Standard error

Example:

```
gene    Name    p       maf     caf     nmiss   ntotal  beta    se      beta.scores     se.scores
ENSG00000096746 10:68337966:T:TG        0.110198283000258       5.09837401075741e-06    5.09837401075741e-06    0       389313  0.739648871041446       0.463061222626501       0.739648871041446       0.463061222626501
```

## Build Jobs script

The `buildJobs.R` script can be used to build a job file for Slurm or other scheduling systems.
This is the best entry point to use the whole pipeline.
The script or its output may need adjustments if you don't use Slurm but other schedulers.

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
SBATCH_ADDITIONAL_PARAMS | optional parameters to pass to 'sbatch' | --exclude=imbi6
MYSQL_EXON_DB | MySQL exon DB, DB name | wuttke
MYSQL_EXON_USER | MySQL exon DB, user name | wuttke
MYSQL_EXON_PASSWORD | MySQL exon DB, password | .....
MYSQL_EXON_HOST | MySQL exon DB, host name | biom131
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

The resulting jobs shell file is organized using Bash functions to group
related jobs. 
You can comment the invocations of these functions at the file end.
This helps to re-run only certain parts of the pipeline, e.g.
if only plotting parameters have been changed.

## BioMart helper files

Two helper files are required for annotation and plotting.

### Genes Mapping File

Use the script `biomart/export_mart_genes.sh`.

This script uses the BioMart web service: http://www.ensembl.org/biomart/martservice

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

Query:
```xml
<Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >
 <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
  <Attribute name = "ensembl_gene_id" />
  <Attribute name = "exon_chrom_start" />
  <Attribute name = "exon_chrom_end" />
  <Attribute name = "ensembl_exon_id" />
  <Attribute name = "rank" />
 </Dataset>
</Query>
```

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

As SQLite is known to have issues with DB files stored on NFS shares,
you can also use a MySQL database.
Helpful scripts are provided in the `biomart` folder.
