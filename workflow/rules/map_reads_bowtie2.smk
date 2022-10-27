rule map_paired_reads:
    input:
        r1 = scratch_dict["trimmed_reads"] / "{sample}.1_trimmed.fastq.gz",
        r2 = scratch_dict["trimmed_reads"] / "{sample}.2_trimmed.fastq.gz",
        index = scratch_dict["done_files"]["index"],
    output:
        temp(scratch_dict["map_reads_paired"] / "{sample}.sam")
    conda:
        "../envs/mapping.yaml"
    log: 
        "logs/mapping/map_paired_reads/{sample}.log"
    # benchmark: 
    #     "benchmark/mapping/map_paired_reads/{sample}.benchmark"
    params:
        name=config["experiment name"],
        index=scratch_dict["genome_index"],
        options=config["read mapping"]["bowtie2 options string"]
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    shell:
        "bowtie2 {params.options} -x {params.index}/{params.name} -1 {input.r1} -2 {input.r2} -p {resources.ntasks} -S {output} &> {log}"

rule map_unpaired_reads:
    input:
        upr = scratch_dict["trimmed_reads"] / "{sample}.unpaired_trimmed.fastq.gz",
        index = scratch_dict["done_files"]["index"],
    output:
        temp(scratch_dict["map_reads_unpaired"] / "{sample}.sam")
    conda:
        "../envs/mapping.yaml"
    log: 
        "logs/mapping/map_unpaired_reads/{sample}.log"
    # benchmark: 
    #     "benchmark/mapping/map_unpaired_reads/{sample}.benchmark"
    params:
        name=config["experiment name"],
        index=scratch_dict["genome_index"],
        options=config["read mapping"]["bowtie2 options string"]
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    shell:
        "bowtie2 {params.options} -x {params.index}/{params.name} -U {input.upr} -p {resources.ntasks} -S {output} &> {log}"

rule sort_paired_bam:
    input:
        scratch_dict["map_reads_paired"] / "{sample}.sam"
    output:
        scratch_dict["map_reads_paired"] / "{sample}.bam"
    conda:
        "../envs/mapping.yaml"
    log: 
        "logs/mapping/sort_paired_reads/{sample}.log"
    benchmark: 
        "benchmark/mapping/sort_paired_reads/{sample}.benchmark"
    params:
        name=config["experiment name"],
        index=scratch_dict["genome_index"],
        options=config["read mapping"]["bowtie2 options string"]
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '64G',
        ntasks = 16,
        time = '12:00:00',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    shell:
        "samtools sort -@ {resources.ntasks} -T {output}.temp -o {output} {input} &> {log}"


rule sort_unpaired_bam:
    input:
        scratch_dict["map_reads_unpaired"] / "{sample}.sam"
    output:
        scratch_dict["map_reads_unpaired"] / "{sample}.bam"
    conda:
        "../envs/mapping.yaml"
    log: 
        "logs/mapping/sort_unpaired_reads/{sample}.log"
    benchmark: 
        "benchmark/mapping/sort_unpaired_reads/{sample}.benchmark"
    params:
        name=config["experiment name"],
        index=scratch_dict["genome_index"],
        options=config["read mapping"]["bowtie2 options string"]
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '64G',
        ntasks = 16,
        time = '12:00:00',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    shell:
        "samtools sort -@ {resources.ntasks} -T {output}.temp -o {output} {input} &> {log}"
