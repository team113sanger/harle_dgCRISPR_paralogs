#!/usr/bin/env bash

# Check repository path exists or exit
if [[ ! -e ${REPO_PATH} ]]
then
	echo "Repository path does not exist: ${REPO_PATH}"
	exit 1;
fi

echo "Repository path: ${REPO_PATH}"

# Function for exiting on error
exit_on_error () {
	err_path=$1
	if [[ -e $err_path ]]
	then
		if [[ ! -z $(grep '[^[:space:]]' $err_path) ]]
		then
			echo "Error found in: ${err_path}"
			exit 1;
		fi
	fi
}

Rscript ${REPO_PATH}/SCRIPTS/dual_guide/Bassik/03_run_bassik_analysis.R \
    -f ${REPO_PATH}/DATA/preprocessing/lfc_matrix.scaled.tsv \
    -y ${REPO_PATH}/DATA/dual_guide/Bassik/pred_vs_obs_y12.tsv \
    -m ${REPO_PATH}/DATA/dual_guide/Bassik/missing_data.tsv \
    -s ${REPO_PATH}/METADATA/sample_annotations.tsv \
    --annotations 13 \
    -c "F1,F2,F3,F4,F5,F6,F7,F8,F9,F10" \
    -t ${REPO_PATH}/DATA/dual_guide/Bassik \
    -r ${REPO_PATH}/DATA/RDS/dual_guide/Bassik \
    2>${REPO_PATH}/LOGS/dual_guide/Bassik/run_bassik_analysis.e \
    1>${REPO_PATH}/LOGS/dual_guide/Bassik/run_bassik_analysis.o
	exit_on_error ${LOG_PATH}/run_bassik_analysis.e
    

echo "DONE."