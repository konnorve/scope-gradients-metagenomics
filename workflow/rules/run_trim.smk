rule run_trim_paired_reads:
    input:
        r1 = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'fwd_read'],
        r2 = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'rev_read'],
        ref = Path(config["input"]["adapter_file"]),
    output:
        o1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq.gz",
        o2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq.gz",
    conda:
        "../envs/bbtools.yaml"
    log:
        "logs/run_trim/run_trim_paired_reads/{sample}.log"
    benchmark:
        "benchmark/run_trim/run_trim_paired_reads/{sample}.benchmark"
    threads: 10
    resources:
        mem_mb = 100000
    shell:
        "bbduk.sh "
        "threads={threads} "
        "in1={input.r1} in2={input.r2} "
        "out1={output.o1} out2={output.o2} "
        "minlen=25 qtrim=rl trimq=10 "
        "ref={input.ref} ktrim=r k=23 mink=11 hdist=1 > {log}"

rule run_trim_unpaired_reads:
    input:
        up = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'unpaired_read'],
        ref = Path(config["input"]["adapter_file"]),
    output:
        scratch_dict["trimmed_reads"] / "{sample}_unpaired_trimmed.fastq.gz",
    conda:
        "../envs/bbtools.yaml"
    log:
        "logs/run_trim/run_trim_unpaired_reads/{sample}.log"
    benchmark:
        "benchmark/run_trim/run_trim_unpaired_reads/{sample}.benchmark"
    threads: 10
    resources:
        mem_mb = 100000
    shell:
        "bbduk.sh "
        "threads={threads} "
        "in={input.up} "
        "out={output} "
        "minlen=25 qtrim=rl trimq=10 "
        "ref={input.ref} ktrim=r k=23 mink=11 hdist=1 > {log}"


