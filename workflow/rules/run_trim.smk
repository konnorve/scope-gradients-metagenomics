rule run_trim_paired_reads:
    input:
        r1 = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'fwd_read'],
        r2 = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'rev_read'],
        ref = Path(config["input"]["adapter file"]),
    output:
        o1 = scratch_dict["trimmed_reads"] / "{sample}.1_trimmed.fastq.gz",
        o2 = scratch_dict["trimmed_reads"] / "{sample}.2_trimmed.fastq.gz",
    conda:
        "../envs/bbtools.yaml"
    log:
        "logs/run_trim/run_trim_paired_reads/{sample}.log"
    # benchmark:
    #     "benchmark/run_trim/run_trim_paired_reads/{sample}.benchmark"
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '64G',
        ntasks = 16,
        time = '12:00:00',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    shell:
        "bbduk.sh "
        "threads={resources.ntasks} "
        "in1={input.r1} in2={input.r2} "
        "out1={output.o1} out2={output.o2} "
        "minlen=25 qtrim=rl trimq=10 "
        "ref={input.ref} ktrim=r k=23 mink=11 hdist=1 &> {log}"

rule run_trim_unpaired_reads:
    input:
        upr = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'unpaired_read'],
        ref = Path(config["input"]["adapter file"]),
    output:
        scratch_dict["trimmed_reads"] / "{sample}.unpaired_trimmed.fastq.gz",
    conda:
        "../envs/bbtools.yaml"
    log:
        "logs/run_trim/run_trim_unpaired_reads/{sample}.log"
    # benchmark:
    #     "benchmark/run_trim/run_trim_unpaired_reads/{sample}.benchmark"
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '64G',
        ntasks = 16,
        time = '12:00:00',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    shell:
        "bbduk.sh "
        "threads={resources.ntasks} "
        "in={input.upr} "
        "out={output} "
        "minlen=25 qtrim=rl trimq=10 "
        "ref={input.ref} ktrim=r k=23 mink=11 hdist=1 &> {log}"
