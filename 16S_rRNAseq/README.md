Pre-processing 16S rRNAseq data using nf-core/taxprofiler and kraken2 database

In Puhti, set working directory: config, DATA, DB, LOGS, RESULTS, SCRIPTS, workflow
- $ cd /scratch/project_XXXX/human_urine

**File upload**
- upload the fastq.qz files into the cd /scratch/project_XXXX/human_urine/DATA folder (e.g. using WinSCP)

**download databases:**
- $ cd  /scratch/project_XXXX/human_urine/DB/

**Human reference database:**
download FASTA file from here: https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_009914755.4/ If you downloaded it from local dir copy it to puhti using scp 
- $ unzip T2T-CHM13v2.0.zip 

**Kraken2**
- $ wget https://genome-idx.s3.amazonaws.com/kraken/16S_Silva138_20200326.tgz 
- $ tar -xzf 16S_Silva138_20200326.tar.gz

 **Create database.csv**
- $ cd /scratch/project_XXXX/human_urine/config
- Cretae the database sheet: "config/database.csv"
        - Create in Excel, save as .CSV, move to Putty. 

**Create samplesheet**
   - Create the samplesheet: "config/samplesheet.csv"

**Create and run Nextflow scripts**
- $ cd /scratch/project_XXXX/human_urine/SCRIPTS
- $ module load nextflow
- $ nextflow pull nf-core/taxprofiler
- $ nano taxprofiler.sh
- $ sbatch taxprofiler.sh

**Results**
  - Results can be found in "kraken2_Kraken2_combined_reports".
  - Copy this to text file to a .csv file (Excel), and keep only the rows with bacteria
  - The metadata will be used in our Metadata file, so make sure you check to which samples S1, S2, S3 etc belong

  - These results can then be used for downstream processing, using the /SCRIPTS/Downstream_processing_Rstudio


