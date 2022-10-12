# map sample reads to concatenated genome
rule map_filtered_interal_standard_reads:
    input:
        r1 = scratch_dict["kaiju"]["merged_reads"] / config["internal_standard_taxa"] / "{sample}_1_unmerged.fastq.gz",
        r2 = scratch_dict["kaiju"]["merged_reads"] / config["internal_standard_taxa"] / "{sample}_2_unmerged.fastq.gz",
        merged = scratch_dict["kaiju"]["merged_reads"] / config["internal_standard_taxa"] / "{sample}_merged.fastq.gz",
        unpaired = scratch_dict["kaiju"]["taxon_assigned_reads"] / config["internal_standard_taxa"] / "{sample}_unpaired_trimmed.fastq.gz",
        ref = Path(config["input"]["internal_standard_ref"]),
        indexing = scratch_dict["done_files"]["index"]
    output:
        sam_out = temp(scratch_dict["read_mapping"] / "filtered_mapping" / "{sample}_mapped.sam"),
    conda:
        "../envs/bowtie2.yaml"
    log: 
        "logs/bowtie2/map_filtered_interal_standard_reads/{sample}.log"
    benchmark: 
        "benchmark/bowtie2/map_filtered_interal_standard_reads/{sample}.benchmark"
    shell:
        "bowtie2 -x {input.ref} -1 {input.r1} -2 {input.r2} -U {input.merged},{input.unpaired} -S {output.sam_out} &> {log}"

rule map_unfiltered_interal_standard_reads:
    input:
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq.gz",
        r2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq.gz",
        unpaired = scratch_dict["trimmed_reads"] / "{sample}_unpaired_trimmed.fastq.gz",
        ref = Path(config["input"]["internal_standard_ref"]),
        indexing = scratch_dict["done_files"]["index"]
    output:
        sam_out = temp(scratch_dict["read_mapping"] / "unfiltered_mapping" / "{sample}_mapped.sam"),
    conda:
        "../envs/bowtie2.yaml"
    log: 
        "logs/bowtie2/map_unfiltered_interal_standard_reads/{sample}.log"
    benchmark: 
        "benchmark/bowtie2/map_unfiltered_interal_standard_reads/{sample}.benchmark"
    shell:
        "bowtie2 -x {input.ref} -1 {input.r1} -2 {input.r2} -U {input.unpaired} -S {output.sam_out} &> {log}"

rule index_genome:
    input:
        ref = Path(config["input"]["internal_standard_ref"])
    output:
        touch(scratch_dict["done_files"]["index"]),
    conda:
        "../envs/bowtie2.yaml"
    log: 
        "logs/bowtie2/index_genome.log"
    benchmark: 
        "benchmark/bowtie2/index_genome.benchmark"
    shell:
        "bowtie2-build {input.ref} {input.ref} &> {log}"

