#BSUB -q normal
#BSUB -J pred_vs_obs[1-800]
#BSUB -oo "/lustre/scratch124/casm/team113/users/vo1/5429_paralog_sl_vicky_finalised_for_paper/LOGS/dual_guide/Bassik/get_pred_vs_obs.%I.o"
#BSUB -eo "/lustre/scratch124/casm/team113/users/vo1/5429_paralog_sl_vicky_finalised_for_paper/LOGS/dual_guide/Bassik/get_pred_vs_obs.%I.e"
#BSUB -R "select[mem>2000] rusage[mem=2000]"
#BSUB -M 2000

# Set repository path
repo_path="/lustre/scratch124/casm/team113/users/vo1/5429_paralog_sl_vicky_finalised_for_paper"

# Set number of jobs (should mirror array range in LSF header)
njobs=800

# Set inputs and output paths
lfc_matrix="${repo_path}/DATA/preprocessing/lfc_matrix.scaled.tsv"
output_directory="${repo_path}/DATA/dual_guide/Bassik"
log_directory="${repo_path}/LOGS/dual_guide/Bassik"

# Load dual guide matrix (dual guides mapped to their corresponding singles )
doubles_gm_file="${repo_path}/METADATA/libraries/dual_guide_matrix.tsv"
doubles_gm_rownum=$(wc -l "${doubles_gm_file}" | awk '{print $1-1}')

# Calculate chunk size from length of dual guide matrix and number of jobs
chunk_size=$(Rscript -e "ceiling( ${doubles_gm_rownum} / ${njobs} ) " | awk '{print $2}' )

# Prepare command to get predicted and observed LFCs for chunk
obs_pred_cmd="Rscript ${repo_path}/SCRIPTS/dual_guide/Bassik/01_get_obs_vs_pred_values_by_chunk.R \
  --fc ${lfc_matrix} \
  --annotations 13 \
  --doubles_guide_matrix ${doubles_gm_file} \
  --helper ${repo_path}/SCRIPTS/dual_guide/helper.R \
  -o ${output_directory} \
  -n ${chunk_size} \
  -i ${LSB_JOBINDEX}"

# Clean out previous results / logs
#if [[ "$(ls -A ${log_directory})" ]]
#then
#	rm ${log_directory}/get_pred_vs_obs.[0-9]*.log
#	rm ${log_directory}/get_pred_vs_obs_y12.[0-9]*.o
#	rm ${log_directory}/get_pred_vs_obs_y12.[0-9]*.e
#	rm ${log_directory}/get_pred_vs_obs_y12.log
#	rm ${log_directory}/get_pred_vs_obs_y12.o
#	rm ${log_directory}/get_pred_vs_obs_y12.e
#	rm ${log_directory}/list_of_pred_vs_obs_y12_files.txt
#	rm ${log_directory}/list_of_missing_data_files.txt
#fi

#if [[ "$(ls -A ${output_directory})" ]]
#then
#	rm ${output_directory}/missing_data.[0-9]*.tsv
#	rm ${output_directory}/pred_vs_obs_y12.[0-9]*.tsv
#	rm ${output_directory}/missing_data.tsv
#	rm ${output_directory}/pred_vs_obs_y12.tsv
#fi 

# Run script
echo "Getting observed and predicted values (${LSB_JOBINDEX} of ${njobs}): "$(date +"%T")
echo "${obs_pred_cmd}"
eval "$obs_pred_cmd > ${log_directory}/get_pred_vs_obs_y12.${LSB_JOBINDEX}.log 2>&1 &"