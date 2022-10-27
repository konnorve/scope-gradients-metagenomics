
rule index_reference_genomes:
    input:
        ref_db=Path(config["input"]["reference genomes"])
    output:
        out_index=directory(scratch_dict["genome_index"]),
        done_file=touch(scratch_dict["done_files"]["index"]),
    conda:
        "../envs/mapping.yaml"
    log: 
        "logs/ref_db_agg/index_genome.log"
    # benchmark: 
    #     "benchmark/ref_db_agg/index_genome.benchmark"
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '1-0',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    params:
        name=config["experiment name"],
        smk_wd=Path('.').absolute()
    shell:
        """
        mkdir {output.out_index}
        cd {input.ref_db}
        bowtie2-build --threads {resources.ntasks} $(ls *.fna | tr '\n' ',') {output.out_index}/{params.name} &> {params.smk_wd}/{log}
        """

rule concat_gff:
    input:
        config["input"]["reference annotations"]
    output:
        scratch_dict["gff_concat"]
    log: 
        "logs/ref_db_agg/concat_gff.log"
    # benchmark: 
    #     "benchmark/ref_db_agg/concat_gff.benchmark"
    resources: 
        partition = 'sched_mit_chisholm',
        mem = '64G',
        ntasks = 16,
        time = '12:00:00',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    shell:
        """
        echo "##gff-version 3" > {output}
        cat {input}/*.gff | grep -v '^#' >> {output}
        """