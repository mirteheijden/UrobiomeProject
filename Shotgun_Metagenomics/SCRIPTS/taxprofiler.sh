#!/bin/bash
#SBATCH --time=60:00:00
#SBATCH --partition=small
#SBATCH --account=project_2008347
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=80G

export SINGULARITY_TMPDIR=$PWD
export SINGULARITY_CACHEDIR=$PWD
unset XDG_RUNTIME_DIR

# Activate Nextflow on Puhti
module load nextflow

# Run the nextflow pipeline
nextflow run nf-core/taxprofiler -r 1.1.5 -resume \
   -profile singularity --max_cpus 4 \
   --input ./config/samplesheet.csv \
   --databases ./config/database.csv \
   --outdir ./RESULTS \
   --perform_shortread_qc --save_preprocessed_reads \
   --perform_shortread_complexityfilter --shortread_complexityfilter_tool bbduk \
   --perform_shortread_hostremoval \
   --hostremoval_reference /scratch/project_2008347/nextflow_metagenomics/DB/human_CHM13/human_CHM13/chm13.draft_v1.0_plusY.fasta \
   --shortread_hostremoval_index /scratch/project_2008347/nextflow_metagenomics/DB/human_CHM13/human_CHM13/ \
   --save_hostremoval_unmapped \
   --save_analysis_ready_fastqs \
   --run_profile_standardisation \
   --run_motus --run_metaphlan
