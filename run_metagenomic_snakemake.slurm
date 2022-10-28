#!/usr/bin/env bash
#SBATCH --job-name=metasmk
#SBATCH --time 10-0                         # wall clock limit
#SBATCH -p sched_mit_chisholm               # partition selection
#SBATCH -n 1                                # tasks (essentially threads)
#SBATCH -o logs/%j_slurm_output.txt
#SBATCH -e logs/%j_slurm_error.txt

source activate snakemake

SNAKEFILE=workflow/Snakefile

# lints and formats snakefile
# snakemake --snakefile $SNAKEFILE --lint

snakemake --snakefile $SNAKEFILE --unlock

# run job:
snakemake \
    --snakefile $SNAKEFILE \
    --use-conda \
    --conda-frontend conda \
    --rerun-incomplete \
    --cores all
    # --cluster "sbatch -p {resources.partition} -t {resources.time} -n {resources.ntasks} --mem={resources.mem} -o {resources.output} -e {resources.error}" \
    # --jobs 1 \
    # --keep-going

# get report once job is finished
# snakemake --snakefile $SNAKEFILE --report report.html