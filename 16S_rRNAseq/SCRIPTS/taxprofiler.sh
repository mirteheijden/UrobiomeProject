#!/bin/bash
#SBATCH --job-name=taxprofiler
#SBATCH --output=taxprofiler.%j.out
#SBATCH --error=taxprofiler.%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --account=project_200XXXX
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G

export SINGULARITY_TMPDIR=$PWD
export SINGULARITY_CACHEDIR=$PWD
unset XDG_RUNTIME_DIR

# Activate  Nextflow on Puhti
module load nextflow
module load minimap2

#run the nextflow
# Go to the folder containing your input files
cd /scratch/project_200XXXX/human_urine

# Run nf-core/taxprofiler
nextflow run nf-core/taxprofiler -r 1.1.5 -resume \
   -profile singularity --max_cpus 4 \
   --container-engine singularity \
   --input /scratch/project_200XXXX/human_urine/config/samplesheet.csv \
   --databases /scratch/project_200XXXX/human_urine/config/database.csv \
   --outdir ./RESULTS/batch2  \
   --perform_shortread_qc \
   --save_preprocessed_reads \
   --shortread_qc_qualityfilter_keeppercent 90 \
   --perform_shortread_complexityfilter  --shortread_complexityfilter_tool bbduk \
   --perform_shortread_hostremoval \
   --hostremoval_reference /scratch/project_200XXXX/nextflow_metagenomics/DB/human_CHM13/human_CHM13/chm13.draft_v1.0_plusY.fasta \
   --longread_hostremoval_index /scratch/project_200XXXX/nextflow_metagenomics/DB/human_CHM13/host_index_chm13 \
   --save_hostremoval_unmapped \
   --save_analysis_ready_fastqs \
   --run_profile_standardisation \
   --run_kraken2 \
   -with-report -with-trace -with-dag
