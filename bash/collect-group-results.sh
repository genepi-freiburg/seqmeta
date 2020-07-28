#!/bin/bash

## adjust these three variables
## run in parent directory of the output directory
PARAMETER_FILE=parameter-file
SUMMARY_FILE=group-summary.txt
SUMMARY_DB=group-summary.sqlite

SAVE_GROUPS=$GROUPS
unset GROUPS
. $PARAMETER_FILE

rm -f $SUMMARY_FILE

for PHENOTYPE in $(echo $PHENOTYPES | sed "s/,/ /g")
do

	for MY_GROUP in $(echo $GROUPS | sed 's/,/ /g')
	do
		for CHR in `seq 1 22`
		do
			GROUP_OUT_FILE="output/group-${MY_GROUP}-${PHENOTYPE}-chr${CHR}.txt"

			if [ ! -f "$GROUP_OUT_FILE" ]
			then
				echo "ERROR: Group output file not found: $GROUP_OUT_FILE"
				exit 3
			fi

			if [ ! -f "$SUMMARY_FILE" ]
			then
				LINE=`head -n1 $GROUP_OUT_FILE`
				echo "group_model	phenotype	$LINE" > $SUMMARY_FILE
			fi

			cat $GROUP_OUT_FILE | tail -n+2 | awk -vgroup=$MY_GROUP -vpheno=$PHENOTYPE '{ print group "\t" pheno "\t" $0 }' >> $SUMMARY_FILE

		done
		echo "Collected: $PHENOTYPE / $MY_GROUP"
	done

done

GROUPS=$SAVE_GROUPS


rm -f $SUMMARY_DB

echo "Create database"
cat | sqlite3 $SUMMARY_DB <<END
CREATE TABLE gene_summary (
  "group_model" TEXT,
  "phenotype" TEXT,
  "gene" TEXT,
  "skat_p" TEXT,
  "skat_Qmeta" TEXT,
  "skat_cmaf" TEXT,
  "skat_nmiss" TEXT,
  "skat_nsnps" TEXT,
  "burden_p" TEXT,
  "burden_beta" TEXT,
  "burden_se" TEXT,
  "burden_cmafTotal" TEXT,
  "burden_cmafUsed" TEXT,
  "burden_nsnpsTotal" TEXT,
  "burden_nsnpsUsed" TEXT,
  "burden_nmiss" TEXT,
  "skat_o_p" TEXT,
  "skat_o_pmin" TEXT,
  "skat_o_rho" TEXT,
  "skat_o_cmaf" TEXT,
  "skat_o_nmiss" TEXT,
  "skat_o_nsnps" TEXT,
  "skat_o_errflag" TEXT
);
END

echo "Load table and create index"
sqlite3 $SUMMARY_DB ".separator \t" ".import '| tail -n+2 $SUMMARY_FILE' gene_summary" "create index gene on gene_summary (gene);"


