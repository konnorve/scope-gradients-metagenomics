from pathlib import Path
import pandas as pd
import numpy as np; np.set_printoptions(suppress=True)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns; sns.set()

def add_col_color_legends(col_color_dict, fig, ax):
    offset = 0
    for attr, info in col_color_dict.items():
        patches = [mpatches.Patch(facecolor=c, label=l) for l, c in info.items()]
        attr_legend = plt.legend(
                    handles=patches, 
                    title=attr,
                    bbox_to_anchor=(0.95 + offset, 1), 
                    bbox_transform=fig.transFigure, 
                    loc='upper right')
        ax.add_artist(attr_legend)
        offset += 0.05

out_dir = Path(snakemake.output[0])
out_dir.mkdir()

gene_attr = list(snakemake.config['omegas heatmap']['columns'].keys())
phys_attr = list(snakemake.config['omegas heatmap']['rows'].keys())

zscores = pd.read_table(snakemake.input['zscores'], index_col='cycog_iid').T
phys_data = pd.read_table(snakemake.input['physiology_data'], index_col=['cruise','sample_id'])
phys_data = phys_data.sort_values(phys_attr)

col_color_dict = {}
for attr, palette in snakemake.config['omegas heatmap']['columns'].items():
    unique_attrs = set()
    for omega in snakemake.config['omegas'].values():
        if attr in omega['attributes'].keys():
            unique_attrs.add(omega['attributes'][attr])
    if len(unique_attrs.intersection(set(['Low','Medium','High']))) == 3:
        unique_attrs = ['Low','Medium','High']
    if 'light' in palette:
        attr_color_scheme = sns.color_palette(palette, len(unique_attrs)+2)[1:-1]
    else:
        attr_color_scheme = sns.color_palette(palette, len(unique_attrs))
    col_color_dict[attr] = dict(zip(unique_attrs, attr_color_scheme))

data_frames = {}
gene_dict = {}
for omega in snakemake.config['omegas'].values():
    for g, cy in omega['genes'].items():
        gene_dict[cy] = g
    data_frames[tuple(map(omega['attributes'].get, gene_attr))] = zscores[list(set(omega['genes'].values()).intersection(set(zscores.columns)))]

zscores = pd.concat(data_frames, axis='columns', names=gene_attr)

col_colors = [list(map(color_dict.get, zscores.columns.get_level_values(attr).to_list())) for attr, color_dict in col_color_dict.items()]

for cruise in phys_data.index.unique(level='cruise'):
    phys_subset = phys_data.loc[cruise, ][phys_attr]
    z_subset = zscores.join(phys_subset, how='inner')
    z_subset = z_subset.set_index(phys_subset.columns.to_list(), append=True)
    z_subset.columns = pd.MultiIndex.from_tuples(z_subset.columns, names=gene_attr+['cycog'])

    row_color_dict = {}
    for attr, palette in snakemake.config['omegas heatmap']['rows'].items():
        unique_attrs = phys_data.loc[cruise][attr].sort_values().to_list()
        if 'light' in palette:
            attr_color_scheme = sns.color_palette(palette, len(unique_attrs)+2)[1:-1]
        else:
            attr_color_scheme = sns.color_palette(palette, len(unique_attrs))
        row_color_dict[attr] = dict(zip(unique_attrs, attr_color_scheme))

    row_colors = [list(map(color_dict.get, z_subset.index.get_level_values(attr).to_list())) for attr, color_dict in row_color_dict.items()]
    row_colors = [[c if c else (1,1,1) for c in clist] for clist in row_colors]

    z_subset.index = z_subset.index.get_level_values('latitude')
    z_subset.columns = list(map(gene_dict.get, z_subset.columns.get_level_values('cycog').to_list()))
    
    fig = sns.clustermap(z_subset, figsize=(30,10), col_cluster=False, row_cluster=False, col_colors=col_colors, row_colors=row_colors, dendrogram_ratio=0.15)
    add_col_color_legends(col_color_dict | row_color_dict, plt.gcf(), plt.gca())
    fig.savefig(out_dir / f"{snakemake.wildcards['organism']}_{cruise}_no_clust.png")
    plt.close()

    fig = sns.clustermap(z_subset, figsize=(30,10), col_cluster=False, row_cluster=True, col_colors=col_colors, row_colors=row_colors, dendrogram_ratio=0.15)
    add_col_color_legends(col_color_dict | row_color_dict, plt.gcf(), plt.gca())
    fig.savefig(out_dir / f"{snakemake.wildcards['organism']}_{cruise}_row_clust.png")
    plt.close()

    fig = sns.clustermap(z_subset, figsize=(30,10), col_cluster=False, row_cluster=True, col_colors=col_colors, dendrogram_ratio=0.15)
    add_col_color_legends(col_color_dict, plt.gcf(), plt.gca())
    fig.savefig(out_dir / f"{snakemake.wildcards['organism']}_{cruise}_row_clust_no_color.png")
    plt.close()

    fig = sns.clustermap(z_subset, figsize=(30,10), col_cluster=True, row_cluster=True, col_colors=col_colors, row_colors=row_colors, dendrogram_ratio=0.15)
    add_col_color_legends(col_color_dict | row_color_dict, plt.gcf(), plt.gca())
    fig.savefig(out_dir / f"{snakemake.wildcards['organism']}_{cruise}_all_clust.png")
    plt.close()
