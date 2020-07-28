#!/bin/bash

SAVE_GROUPS=$GROUPS
unset GROUPS
. parameter-file

for PHENOTYPE in $(echo $PHENOTYPES | sed "s/,/ /g")
do

	for MY_GROUP in $(echo $GROUPS | sed 's/,/ /g')
	do
		for CHR in `seq 1 22`
		do
			echo -n "Check: $PHENOTYPE, $MY_GROUP, chr$CHR: "

			LOG_FILE="logs/${MY_GROUP}_${PHENOTYPE}-chr${CHR}-*.log"
			LOG_FILE_RESOLVED=`ls $LOG_FILE`
			if [ ! -f "$LOG_FILE_RESOLVED" ]
			then
				echo "ERROR: Log file not found: $LOG_FILE"
				exit 1
			else
				OK=`grep Finished $LOG_FILE_RESOLVED`
				if [ "$OK" == "" ]
				then
					echo "Log file looks incomplete: $LOG_FILE_RESOLVED"
					exit 2
				fi
			fi

			# group-ultraRare-egfr_ckdepi_creat_cys-chr2.txt  sv-syn-bun_serum-chr9.txt
			GROUP_OUT_FILE="output/group-${MY_GROUP}-${PHENOTYPE}-chr${CHR}.txt"
			SV_OUT_FILE="output/sv-${MY_GROUP}-${PHENOTYPE}-chr${CHR}.txt"

			if [ ! -f "$GROUP_OUT_FILE" ]
			then
				echo "ERROR: Group output file not found: $GROUP_OUT_FILE"
				exit 3
			fi

			if [ ! -f "$SV_OUT_FILE" ]
			then
				echo "ERROR: Single-variant output file not found: $SV_OUT_FILE"
				exit 4
			fi

			echo "OK."
		done

		GENE_PLOT_FILE_PATTERN="output/genePlot_${MY_GROUP}-${PHENOTYPE}_*.pdf"
		GENE_PLOT_FILES=`ls $GENE_PLOT_FILE_PATTERN`
		if [ "$GENE_PLOT_FILES" == "" ]
		then
			echo "WARNING: No gene plots found for ${MY_GROUP} and ${PHENOTYPE}!"
		fi

		QQ_PLOT="output/$MY_GROUP-$PHENOTYPE-qq.pdf"
		if [ ! -f "$QQ_PLOT" ]
		then
			echo "ERROR: QQ plot missing for ${MY_GROUP} and ${PHENOTYPE}: ${QQ_PLOT}!"
			exit 5
                fi

                TOP_FILE="output/$MY_GROUP-$PHENOTYPE-top.txt"
                if [ ! -f "$TOP_FILE" ]
                then
                        echo "ERROR: Top genes file missing for ${MY_GROUP} and ${PHENOTYPE}: ${TOP_FILE}!"
			exit 6
                fi

                XLSX_FILE="output/$MY_GROUP-$PHENOTYPE-top.xlsx"
                if [ ! -f "$XLSX_FILE" ]
                then
                        echo "ERROR: XLSX file missing for ${MY_GROUP} and ${PHENOTYPE}: ${XLSX_FILE}!"
			exit 7
                fi
	done

done

GROUPS=$SAVE_GROUPS
