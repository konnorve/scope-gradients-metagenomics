#!/usr/bin/env bash
#SBATCH --job-name=metasmk
#SBATCH --time 2-0                         # wall clock limit
#SBATCH -p sched_mit_chisholm               # partition selection
#SBATCH -n 1                                # tasks (essentially threads)
#SBATCH -o logs/%j_slurm_output.txt
#SBATCH -e logs/%j_slurm_error.txt

source activate snakemake

# lints and formats snakefile
# snakemake --snakefile workflow/Snakefile --lint

snakemake --snakefile workflow/Snakefile --unlock

# run job:
snakemake \
    --snakefile workflow/Snakefile \
    --use-conda \
    --conda-frontend conda \
    --rerun-incomplete \
    --cluster "sbatch -p {resources.partition} -t {resources.time} -n {resources.ntasks} --mem={resources.mem} -o {resources.output} -e {resources.error}" \
    --jobs 5 \
    --keep-going

# get report once job is finished
# snakemake --snakefile $SNAKEFILE --report report.html