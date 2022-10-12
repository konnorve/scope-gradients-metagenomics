#!/usr/bin/env bash
# 
#SBATCH --job-name=metasmk
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 10-0
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
#SBATCH --comment="Metagenome Analysis Workflow"
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
# #SBATCH -N 2
#SBATCH -n 40
# #SBATCH --mem 250G
#
#SBATCH -o logs/%j_slurm_output.txt
#SBATCH -e logs/%j_slurm_error.txt

source activate snakemake

SNAKEFILE=workflow/Snakefile

# lints and formats snakefile
snakemake --snakefile $SNAKEFILE --lint

snakemake --snakefile $SNAKEFILE --unlock

# run job:
snakemake --cores all --resources mem_mb=500000 --snakefile $SNAKEFILE --use-conda --conda-frontend conda --rerun-incomplete

# get report once job is finished
# snakemake --snakefile $SNAKEFILE --report report.html