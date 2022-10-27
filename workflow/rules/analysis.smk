

rule generate_org_count_tables:
    input:
        htseq_out=expand(scratch_dict["counting_dir"] / "{sample}.tsv", sample=SAMPLES),
        sccg_list=lambda wildcards: config['input']['sccg lists'][wildcards.organism],
        img_lookup_table=config['input']['reference genome lookup table'],
        cycog_lookup_table=config['input']['cycog lookup table'],
    output:
        results_dict['organism_tables'] / "{organism}_count_table.tsv"
    conda:
        "../envs/analysis.yaml"
    resources:
        partition = 'sched_mit_chisholm',
        mem = '10G',
        ntasks = 1,
        time = '0-1',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    script:
        "../scripts/generate_org_count_table.py"


rule make_gene_pca:
    input:
        organism_count_table = results_dict['organism_tables'] / "{organism}_count_table.tsv",
        nutrient_stress_genes = config['input']['nutrient stress genes'],
        physiology_data = config['input']['physiology data'],
    output:
        zscores = results_dict['zscores'] / "{organism}_zscores.tsv",
        omegas = results_dict['omegas'] / "{organism}_omegas.tsv",
        loadings_df = results_dict['pca_loadings'] / "{organism}_pca_loadings.tsv",
        loadings_fig = results_dict['pca_loadings'] / "{organism}_pca_loadings.html",
        transformed_df = results_dict['pca_overlays'] / "{organism}" / "{organism}_ordination.tsv",
    conda:
        "../envs/analysis.yaml"
    resources:
        partition = 'sched_mit_chisholm',
        mem = '10G',
        ntasks = 1,
        time = '0-1',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    script:
        "../scripts/make_gene_pca.py"


rule plot_pca_overlays:
    input:
        results_dict['pca_overlays'] / "{organism}" / "{organism}_ordination.tsv"
    output:
        results_dict['pca_overlays'] / "{organism}" / "{organism}_{column}_pca_overlay.png"
    conda:
        "../envs/analysis.yaml"
    resources:
        partition = 'sched_mit_chisholm',
        mem = '10G',
        ntasks = 1,
        time = '0-1',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    script:
        "../scripts/plot_pca_overlays.py"

rule make_omega_regression:
    input:
        omegas = results_dict['omegas'] / "{organism}_omegas.tsv",
        physiology_data = config['input']['physiology data'],
    output:
        directory(results_dict['regressions'] / "{organism}")
    conda:
        "../envs/analysis.yaml"
    resources:
        partition = 'sched_mit_chisholm',
        mem = '10G',
        ntasks = 1,
        time = '0-1',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    script:
        "../scripts/make_omega_regression.py"

rule make_gene_heatmap:
    input:
        zscores = results_dict['zscores'] / "{organism}_zscores.tsv",
        nutrient_stress_genes = config['input']['nutrient stress genes'],
        physiology_data = config['input']['physiology data'],
    output:
        directory(results_dict['heatmaps'] / "{organism}")
    conda:
        "../envs/analysis.yaml"
    resources:
        partition = 'sched_mit_chisholm',
        mem = '10G',
        ntasks = 1,
        time = '0-1',
        output = 'logs/smk_slurm/%j_slurm.out',
        error = 'logs/smk_slurm/%j_slurm.err',
    script:
        "../scripts/make_gene_heatmap.py"

