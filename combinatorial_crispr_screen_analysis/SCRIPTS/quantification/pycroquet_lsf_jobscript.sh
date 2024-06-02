#BSUB -q long
#BSUB -J pycroquet[1-516]
#BSUB -oo "LOGS/pyCROQUET/pycroquet.%I.o"
#BSUB -eo "LOGS/pyCROQUET/pycroquet.%I.e"
#BSUB -R "select[mem>10000] rusage[mem=10000] span[hosts=1]"
#BSUB -M 10000
#BSUB -n 8

module load pycroquet/1.5.1

guides="${REPO_PATH}/METADATA/libraries/paralog_library.pycroquet.tsv"
cram_directory="${REPO_PATH}/DATA/sequencing"
cram_files=($(ls ${cram_directory}/*.cram))
array_index=$(expr ${LSB_JOBINDEX} - 1)
query_file_name=$(basename ${cram_files[${array_index}]})
query_file_label=${query_file_name::-5}
output_directory="${REPO_PATH}/DATA/pyCROQUET"
 
pycroquet dual-guide -g ${guides} -q "${cram_directory}/${query_file_name}" -o "${output_directory}/${query_file_label}" --chunks 50000 -b exact -c 8
