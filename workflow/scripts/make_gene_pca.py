import pandas as pd
from sklearn.decomposition import PCA
import plotly.graph_objects as go
import plotly.express as px
import numpy as np

counts = pd.read_table(snakemake.input['organism_count_table'], index_col="cycog_iid")
nutrient_concentrations = pd.read_table(snakemake.input['physiology_data'], index_col="sample_id")

cycogs_in_config = set()
for v in snakemake.config['omegas'].values():
    cycogs_in_config |= set(v['genes'].values())
cycogs_in_counts = set(counts.index.to_list())
cycogs_for_analysis = cycogs_in_config.intersection(cycogs_in_counts)

counts = counts.loc[list(cycogs_for_analysis)]

# convert frequencies into z scores
zscores = counts.sub(counts.mean(axis=1), axis=0).div(counts.std(axis=1), axis=0)
zscores.to_csv(snakemake.output['zscores'], sep='\t')

omegas = []
for omega, info in snakemake.config['omegas'].items():
    omega_series = zscores.loc[list(set(info['genes'].values()).intersection(cycogs_for_analysis))].mean()
    omega_series.name = omega
    omegas.append(omega_series)
omegas = pd.concat(omegas, axis=1)
omegas.index.name = "sample_id"
omegas.to_csv(snakemake.output['omegas'], sep='\t')

pca = PCA(n_components=2).fit(zscores.T)

transformed_df = pd.DataFrame(pca.transform(zscores.T), columns=['x', 'y'], index=zscores.columns)
transformed_df = transformed_df.join(nutrient_concentrations)
transformed_df.to_csv(snakemake.output['transformed_df'], sep='\t')

loadings_df = pd.DataFrame(pca.components_.T*np.sqrt(pca.explained_variance_), columns=['x', 'y'], index=zscores.index)
loadings_df.to_csv(snakemake.output['loadings_df'], sep='\t')

fig = go.Figure()

unique_attributes = dict()
for omega, info in snakemake.config['omegas'].items():
    for k, v in info['attributes'].items():
        if k in unique_attributes.keys():
            unique_attributes[k].add(v)
        else:
            unique_attributes[k] = set([v])

# make attribute marker and color dictionary
color_attr = 'nutrient'
marker_attr = 'level'

color_dict = {t:c for t, c in zip(unique_attributes[color_attr], px.colors.qualitative.Plotly)}
marker_dict = {t:c for t, c in zip(unique_attributes[marker_attr], ['circle', 'square', 'diamond', 'cross', 'x'])}

for omega, info in snakemake.config['omegas'].items():
    for gene, cycog in info['genes'].items():
        if cycog in cycogs_for_analysis:
            fig.add_trace(
                go.Scatter(
                    x=[0, loadings_df.loc[cycog, 'x']],
                    y=[0, loadings_df.loc[cycog, 'y']],
                    name=gene,
                    legendgroup=omega,
                    legendgrouptitle_text=omega,
                    marker_color=color_dict[info['attributes'][color_attr]],
                    marker_symbol=marker_dict[info['attributes'][marker_attr]],
                    # text=f"cycog: {cycog}<br>"+"<br>".join([f"{k}: {v}" for k,v in row.to_dict().items()])
                )
            )

for omega, info in snakemake.config['omegas'].items():
    omega_series = loadings_df.loc[list(set(info['genes'].values()).intersection(cycogs_for_analysis))].mean()
    fig.add_trace(
        go.Scatter(
            x=[omega_series.loc['x']],
            y=[omega_series.loc['y']],
            name=f'&#937; {omega}',
            legendgroup="Omegas",
            legendgrouptitle_text="Omegas",
            marker_color='black',
            marker_symbol='circle',
            mode="text+markers",
            text=f'&#937; {omega}',
            textposition="top center"
        )
    )

fig.update_layout(
    xaxis_title=f"PCA1 ({pca.explained_variance_ratio_[0]:.0%} of variance explained)",
    yaxis_title=f"PCA2 ({pca.explained_variance_ratio_[1]:.0%} of variance explained)",
)
fig.write_html(snakemake.output['loadings_fig'])
