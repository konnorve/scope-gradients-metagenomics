
rule count_reads:
    input:
        pe_bam = scratch_dict["map_reads_paired"] / "{sample}.bam",
        se_bam = scratch_dict["map_reads_unpaired"] / "{sample}.bam",
        gff = scratch_dict["gff_concat"]
    output:
        scratch_dict["counting_dir"] / "{sample}.tsv",
    conda:
        "../envs/htseq-count.yaml"
    log: 
        "logs/count_reads/{sample}.log"
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '10G',
        ntasks = 1,
        time = '12:00:00',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    params:
        mapq_thresh=config["counting features"]["minimum mapq quality"]
    shell:
        "htseq-count --idattr=gene -a {params.mapq_thresh} -s no -t CDS -r pos --nonunique all -c {output} {input.pe_bam} {input.se_bam} {input.gff}"