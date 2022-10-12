
rule assign_paired_reads:
    input:
        cycog_db_index = Path(config["input"]["cycog_db_index"]),
        cycog_db_nodes = Path(config["input"]["cycog_db_nodes"]),
        R1 = scratch_dict["kaiju"]["merged_reads"] / "{taxon}" / "{sample}_1_unmerged.fastq.gz",
        R2 = scratch_dict["kaiju"]["merged_reads"] / "{taxon}" / "{sample}_2_unmerged.fastq.gz",
    output:
        scratch_dict["cycog"]["kaiju_assignment"] / "{taxon}" / "{sample}_paired_kaiju.out"
    conda:
        "../envs/kaiju.yaml"
    log:
        "logs/cycog/assign_paired_reads/{taxon}.{sample}.log"
    benchmark:
        "benchmark/cycog/assign_paired_reads/{taxon}.{sample}.benchmark"
    threads: 10
    resources:
        mem_mb = 100000
    shell:
        """
        kaiju -z {threads} \
            -t {input.cycog_db_nodes} \
            -f {input.cycog_db_index} \
            -i {input.R1} \
            -j {input.R2} \
            -o {output} 2> {log}
        """


rule assign_unpaired_reads:
    input:
        cycog_db_index = Path(config["input"]["cycog_db_index"]),
        cycog_db_nodes = Path(config["input"]["cycog_db_nodes"]),
        R = scratch_dict["kaiju"]["merged_reads"] / "{taxon}" / "{sample}_unpaired.fastq.gz",
    output:
        scratch_dict["cycog"]["kaiju_assignment"] / "{taxon}" / "{sample}_unpaired_kaiju.out"
    conda:
        "../envs/kaiju.yaml"
    log:
        "logs/cycog/assign_unpaired_reads/{taxon}.{sample}.log"
    benchmark:
        "benchmark/cycog/assign_unpaired_reads/{taxon}.{sample}.benchmark"
    threads: 10
    resources:
        mem_mb = 100000
    shell:
        """
        kaiju -z {threads} \
            -t {input.cycog_db_nodes} \
            -f {input.cycog_db_index} \
            -i {input.R} \
            -o {output} 2> {log}
        """

rule assign_merged_reads:
    input:
        cycog_db_index = Path(config["input"]["cycog_db_index"]),
        cycog_db_nodes = Path(config["input"]["cycog_db_nodes"]),
        R = scratch_dict["kaiju"]["merged_reads"] / "{taxon}" / "{sample}_merged.fastq.gz",
    output:
        scratch_dict["cycog"]["kaiju_assignment"] / "{taxon}" / "{sample}_merged_kaiju.out"
    conda:
        "../envs/kaiju.yaml"
    log:
        "logs/cycog/assign_merged_reads/{sample}.log"
    benchmark:
        "benchmark/cycog/assign_merged_reads/{sample}.benchmark"
    threads: 10
    resources:
        mem_mb = 100000
    shell:
        """
        kaiju -z {threads} \
            -t {input.cycog_db_nodes} \
            -f {input.cycog_db_index} \
            -i {input.R} \
            -o {output} 2> {log}
        """

rule filter_and_list_taxon_assigned:
    input:
        scratch_dict["cycog"]["kaiju_assignment"] / "{taxon}" / "{sample}_{readpairing}_kaiju.out"
    output:
        scratch_dict["cycog"]["filtering_assignments"] / "{taxon}" / "{gene}" / "{sample}_{readpairing}_kaiju.lst"
    log:
        "logs/cycog/filter_and_list_taxon_assigned/{taxon}.{gene}.{sample}.{readpairing}.log"
    benchmark:
        "benchmark/cycog/filter_and_list_taxon_assigned/{taxon}.{gene}.{sample}.{readpairing}.benchmark"
    params:
        filtr=lambda wildcards, output: config["CyCOG_assignments"][wildcards.gene]
    shell:
        """
        echo "{params.filtr}" > {log}; 
        cat {input} | grep "{params.filtr}" | awk '{{print $2}}' > {output} 2> {log}
        """

rule extract_reads:
    input:
        assignment = scratch_dict["cycog"]["filtering_assignments"] / "{taxon}" / "{gene}" / "{sample}_{readpairing}_kaiju.lst"
        read_group = scratch_dict["kaiju"]["merged_reads"] / "{taxon}" / "{sample}_{readpairing}.fastq.gz",
    output:
        scratch_dict["cycog"]["taxon_assigned_reads"] / "{taxon}" / "{gene}" / "{sample}_{readpairing}.fastq",
    log:
        "logs/cycog/extract_reads/{taxon}.{gene}.{sample}.{readpairing}.log"
    benchmark:
        "benchmark/cycog/extract_reads/{taxon}.{gene}.{sample}.{readpairing}.benchmark"
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk subseq {input.read_group} {input.assignment} > {output} 2> {log}"


rule zip:
    input:
        scratch_dict["cycog"]["taxon_assigned_reads"] / "{taxon}" / "{gene}" / "{anything}.fastq",
    output:
        scratch_dict["cycog"]["taxon_assigned_reads"] / "{taxon}" / "{gene}" / "{anything}.fastq.gz",
    resources:
        mem_mb=100000,
    conda:
        "../envs/gzip.yaml"
    shell:
        "gzip {input}"
