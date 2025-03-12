#!/bin/tcsh
#SBATCH --job-name=Bioseq_Singularity_Test                   # job name
#SBATCH --partition=256GBv2                   # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                             # number of nodes requested by user
#SBATCH --time=02-00:00:00                                    # run time, format: D-H:m:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                             # standard output file name
#SBATCH --error=serialJob.%j.time                             # standard error output file name

### Load the singularity module, 3.9.9 is tested and functional
module load singularity/3.9.9
echo “Singularity loaded”
### This is to use singularity to mount the container image -> run the python script in from the specified directory within the image -> point to relevant input files

while IFS=$'\t' read -r file_name file_loc _; do
    echo "$file_name $file_loc"
    singularity exec /project/GCRB/Hon_lab/shared/container_images/long_read_sequencing/bioseq_v0.sif \
    python /project/GCRB/Hon_lab/s223695/Data_project/20250304_CP3_QC/lrguidetools_3.py \
    -i ${file_loc} \
    -r /project/GCRB/Hon_lab/s223695/Data_project/TFperturb_edist_additonal/dropout_analysis/Pool.total.oligos \
    -n ${file_name}

done < gRNA_seq_information.txt