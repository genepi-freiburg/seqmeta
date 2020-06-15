
SCRIPT_PATH="R/groupTest.R"

mkdir -p output logs

for CHR in `seq 1 22`
do

sbatch $SCRIPT_PATH --chr=$CHR \
	--bgen_path="..../ukb_wes_efe-chr%CHR%.bgen" \
 	--group_file=".../03_rareDamaging/group_file.txt" \
 	--phenotype_file="phenotype.txt" \
 	--phenotype_col="uric_acid_serum" \
 	--covariate_cols="age_crea_serum, sex_male, U1, U2, U3" \
 	--sv_output_path="output/%PHENO%-sv-chr%CHR%.txt" \
 	--group_output_path="output/%PHENO%-group-chr%CHR%.txt" \
	| tee logs/urate-chr$CHR.log

done
