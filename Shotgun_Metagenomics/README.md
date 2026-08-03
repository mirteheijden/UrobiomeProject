Pre-processing shallow shotgun metagenomics sequencing data using nf-core/taxprofiler and MetaPhlan4 database

In Puhti, set working directory: config, DATA, DB, LOGS, RESULTS, SCRIPTS, workflow
- cd /scratch/project_XXXX/human_urine_shotgun_metagenomics

**File upload**
- upload the fastq.qz files into the cd /scratch/project_XXXX/human_urine_shotgun_metagenomics/DATA folder (e.g. using WinSCP)

**download databases:**
- cd  /scratch/project_XXXX/human_urine_shotgun_metagenomics/DB/

**Human reference database:**
download FASTA file from here: https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_009914755.4/ If you downloaded it from local dir copy it to puhti using scp 
- unzip T2T-CHM13v2.0.zip 

**MetaPhlan4 database**
- Follow set-up here: https://docs.csc.fi/apps/metaphlan/ 

 **Create database.csv**
- cd /scratch/project_XXXX/human_urine_shotgun_metagenomics/config
- Cretae the database sheet: "config/database.csv"
        - Create in Excel, save as .CSV, move to Putty. 

**Create samplesheet**
   - Create the samplesheet: "config/samplesheet.csv"

**Create and run Nextflow scripts**
- cd /scratch/project_XXXX/human_urine_shotgun_metagenomics/SCRIPTS
- module load nextflow
- nextflow pull nf-core/taxprofiler
- nano taxprofiler.sh
- sbatch taxprofiler.sh

**Results**
- MetaPhlan4 results "metaphlan_db_meta4_combined_reports.txt" can be found RESULTS/Metaphlan. Use WinSCP to transfer to harddrive or to Rstudio. If they are not merged yet, merge the  .metaphlan_profile.txt files yourself:
- cd /scratch/project_200XXXX/human_feces_batch1/RESULTS/metaphlan
- merge_metaphlan_tables.py db_meta4/*metaphlan_profile.txt > metaphlan_db_meta4_combined_reports.txt
- Check the removal of the human reads by calculating the reads in the Raw data (reads_count_preprocess), after human removal:
- cd ./SCRIPTS
- nano reads_count_preprocess
- nano reads_count_postprocess
- Run the scripts, one by one:
- chmod +x reads_count_preprocess
- ./reads_count_preprocess
- chmod +x reads_count_postprocess
-  ./reads_count_postprocess
-  Results can be found in RESULTS/reads_count_preprocess.tsv, and RESULTS/reads_count_postprocess
-  Move .tsv files to PC using WinSCP


