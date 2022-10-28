import pandas as pd
from pathlib import Path

img_genome_df = pd.read_table(snakemake.input['img_lookup_table'])
cycog_genes = pd.read_table(snakemake.input['cycog_lookup_table'], index_col='gene_iid')

sccg_cycog_list = Path(snakemake.input['sccg_list']).read_text().split("\n")

htseq_metadata_list = []
htseq_count_list = []
for f in snakemake.input['htseq_out']:
    htseq_df = pd.read_table(f, header=0, names=['gene_id', 'paired_counts', 'unpaired_counts'])
    htseq_df.loc[:,'sample_id'] = Path(f).stem

    htseq_counts_df = htseq_df[~htseq_df['gene_id'].str.contains("__")].copy()
    htseq_counts_df.loc[:,'counts'] = htseq_counts_df['paired_counts'] + htseq_counts_df['unpaired_counts']
    htseq_counts_df = htseq_counts_df.join(cycog_genes, on='gene_id')
    htseq_counts_df = htseq_counts_df.set_index(['sample_id','cycog_iid','gene_id'])

    htseq_metadata_df = htseq_df[htseq_df['gene_id'].str.contains("__")].copy()
    htseq_metadata_df.loc[:,'gene_id'] = htseq_metadata_df['gene_id'].str.replace("__", "")
    htseq_metadata_df = htseq_metadata_df.set_index(['sample_id','gene_id'])
    htseq_metadata_df.loc[(Path(f).stem, 'read_counts'), :] = htseq_counts_df.sum(axis=0)

    htseq_metadata_list.append(htseq_metadata_df)
    htseq_count_list.append(htseq_counts_df)

htseq_metadata_df = pd.concat(htseq_metadata_list)
htseq_counts_df = pd.concat(htseq_count_list)

htseq_counts_df = htseq_counts_df[htseq_counts_df['organism']==snakemake.wildcards['organism']]

htseq_counts_df = htseq_counts_df.groupby(level=['sample_id','cycog_iid'])['counts'].sum().unstack(level=0)
sccg_avg_series = htseq_counts_df.loc[sccg_cycog_list].mean(axis=0)
htseq_counts_df = htseq_counts_df.divide(sccg_avg_series, axis='columns')

htseq_counts_df.to_csv(snakemake.output[0], sep='\t')