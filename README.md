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

