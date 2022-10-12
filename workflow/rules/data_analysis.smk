rule plot_bam_depths:
    input:
        scratch_dict["read_mapping"] / "{filter}" / "{sample}_depth.tsv"
    output:
        depth_histogram = results_dict["depth_histograms"] / "{filter}" / "{sample}_depth_histogram.png",
        depth_genome = results_dict["depth_genome"] / "{filter}" / "{sample}_depth_genome.png",
    conda:
        "../envs/deep_variant_calling.yaml"
    script:
        "../scripts/depth_histogram.py"

rule get_fastq_statistics_internal_standards:
    input:
        R1 = expand(scratch_dict["genome_filtered_internal_standards"] / "{{filter}}" / "{sample}_1_unmerged.fastq.gz", sample=SAMPLES),
        R2 = expand(scratch_dict["genome_filtered_internal_standards"] / "{{filter}}" / "{sample}_2_unmerged.fastq.gz", sample=SAMPLES),
        unpaired = expand(scratch_dict["genome_filtered_internal_standards"] / "{{filter}}" / "{sample}_unpaired.fastq.gz", sample=SAMPLES),
    output:
        scratch_dict["genome_filtered_internal_standards"] / "{filter}_read_stats.tsv"
    conda:
        "../envs/seqkit.yaml"
    params:
        ext = ".fastq.gz"
    shell:
        "seqkit stats $(dirname {input.unpaired[0]})/*{params.ext} -T > {output}"

rule get_fastq_statistics_taxons:
    input:
        merged = expand(scratch_dict["kaiju"]["merged_reads"] / "{{taxon}}" / "{sample}_merged.fastq.gz", sample=SAMPLES),
        unmerged_R1 = expand(scratch_dict["kaiju"]["merged_reads"] / "{{taxon}}" / "{sample}_1_unmerged.fastq.gz", sample=SAMPLES),
        unmerged_R2 = expand(scratch_dict["kaiju"]["merged_reads"] / "{{taxon}}" / "{sample}_2_unmerged.fastq.gz", sample=SAMPLES),
        unpaired = expand(scratch_dict["kaiju"]["merged_reads"] / "{{taxon}}" / "{sample}_unpaired.fastq.gz", sample=SAMPLES),
    output:
        scratch_dict["read_stats"] / "{taxon}_read_stats.tsv"
    conda:
        "../envs/seqkit.yaml"
    params:
        ext = ".fastq.gz"
    shell:
        "seqkit stats $(dirname {input.merged[0]})/*{params.ext} -T > {output}"


rule get_fastq_statistics_cycogs:
    input:
        expand(scratch_dict["cycog"]["taxon_assigned_reads"] / "{{taxon}}" / "{{gene}}" / "{sample}_{readpairing}.fastq", sample=SAMPLES, readpairing=["1_unmerged", "2_unmerged", "unpaired", "merged"])
    output:
        scratch_dict["cycog_read_stats"] / "{taxon}" / "{gene}_read_stats.tsv"
    conda:
        "../envs/seqkit.yaml"
    params:
        ext = ".fastq.gz"
    shell:
        "seqkit stats $(dirname {input[0]})/*{params.ext} -T > {output}"


rule make_table:
    input:
        taxon_read_stat_files=expand(scratch_dict["read_stats"] / "{taxon}_read_stats.tsv", taxon=list(config["taxons_of_interst"].keys()))
    output:
        read_count = results_dict['taxon_read_counts'],
        base_count = results_dict['taxon_base_counts']
    conda:
        "../envs/deep_variant_calling.yaml"
    script:
        "../scripts/taxon_read_count_table.py"

rule make_cycog_table:
    input:
        read_stat_files=expand(scratch_dict["cycog_read_stats"] / "{{taxon}}" / "{gene}_read_stats.tsv", gene=list(config["CyCOG_assignments"].keys()))
    output:
        read_count = results_dict['cycog_read_counts'] / "{taxon}_cycog_read_counts.tsv",
        base_count = results_dict['cycog_base_counts'] / "{taxon}_cycog_base_counts.tsv",
    conda:
        "../envs/deep_variant_calling.yaml"
    script:
        "../scripts/taxon_read_count_table.py"

rule collect_internal_standards_stats:
    input:
        only_kaiju=scratch_dict["read_stats"] / f"{config['internal_standard_taxa']}_read_stats.tsv",
        only_bowtie2=scratch_dict["genome_filtered_internal_standards"] / "unfiltered_mapping_read_stats.tsv",
        both_kaiju_bowtie2=scratch_dict["genome_filtered_internal_standards"] / "filtered_mapping_read_stats.tsv"
    output:
        read_count = results_dict['internal_standards_read_counts'],
        base_count = results_dict['internal_standards_base_counts']
    conda:
        "../envs/deep_variant_calling.yaml"
    script:
        "../scripts/internal_standard_counts.py"