
rule convert_sam2bam:
    input:
        scratch_dict["read_mapping"] / "{filter}" / "{sample}_mapped.sam",
    output:
        temp(scratch_dict["read_mapping"] / "{filter}" / "{sample}_mapped.bam"),
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/convert_sam2bam/{sample}.{filter}.log"
    shell:
        "samtools view -S -b {input} > {output} 2> {log}"

rule sort_bam:
    input:
        scratch_dict["read_mapping"] / "{filter}" / "{sample}_mapped.bam",
    output:
        scratch_dict["read_mapping"] / "{filter}" / "{sample}_mapped_sorted.bam"
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/sort_bam/{sample}.{filter}.log"
    shell:
        "samtools sort {input} -o {output} > {log}"

rule bam_coverage:
    input:
        scratch_dict["read_mapping"] / "{filter}" / "{sample}_mapped_sorted.bam"
    output:
        scratch_dict["read_mapping"] / "{filter}" / "{sample}_depth.tsv"
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/samtools_depth/{sample}.{filter}.log"
    shell:
        "samtools depth -o {output} {input}"

rule bam_to_fastq:
    input:
        scratch_dict["read_mapping"] / "{filter}" / "{sample}_mapped_sorted.bam"
    output:
        R1 = scratch_dict["genome_filtered_internal_standards"] / "{filter}" / "{sample}_1_unmerged.fastq.gz",
        R2 = scratch_dict["genome_filtered_internal_standards"] / "{filter}" / "{sample}_2_unmerged.fastq.gz",
        unpaired = scratch_dict["genome_filtered_internal_standards"] / "{filter}" / "{sample}_unpaired.fastq.gz",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/bam_to_fastq/{sample}.{filter}.log"
    threads: 5
    shell:
        "samtools fastq -F 4 -@ {threads} -1 {output.R1} -2 {output.R2} -s {output.unpaired} {input} > {log}"



# rule index_bam:
#     input:
#         scratch_dict["read_mapping"] / "{filter}" / "{sample}_mapped_sorted.bam",
#     output:
#         temp(scratch_dict["read_mapping"] / "{filter}" / "{sample}_mapped_sorted.bam.bai"),
#     resources:
#         mem_mb=100000,
#     conda:
#         "../envs/samtools.yaml"
#     log:
#         "logs/samtools/index_bam/{sample}.log"
#     shell:
#         "samtools index -b {input} > {log}"


# rule index_fasta:
#     input:
#         reference_genome_file,
#     output:
#         reference_genome_index_file,
#     resources:
#         mem_mb=100000,
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         "samtools faidx {input}"
    