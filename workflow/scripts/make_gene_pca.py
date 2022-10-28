import pandas as pd
from sklearn.decomposition import PCA
import plotly.graph_objects as go
import plotly.express as px
import numpy as np

df = pd.read_table(snakemake.input['organism_count_table'], index_col="cycog_iid")
lookup_table = pd.read_table(snakemake.input['nutrient_stress_genes'], index_col="cycog_id")
nutrient_concentrations = pd.read_table(snakemake.input['physiology_data'], index_col="sample_id")

df = df.loc[[cycog for cycog in lookup_table.index.to_list() if cycog in df.index.to_list()]]

# convert frequencies into z scores
df = df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1), axis=0)

df.to_csv(snakemake.output['zscores'], sep='\t')

omegas = df.join(lookup_table[['stress_type', 'stress_level']]).groupby(['stress_type', 'stress_level']).mean()
omegas = omegas.T
omegas.columns = ["_".join(col) for col in omegas.columns.values]
omegas.index.name = "sample_id"
omegas.to_csv(snakemake.output['omegas'], sep='\t')

pca = PCA(n_components=2).fit(df.T)

transformed_df = pd.DataFrame(pca.transform(df.T), columns=['x', 'y'], index=df.columns)
transformed_df = transformed_df.join(nutrient_concentrations)
transformed_df.to_csv(snakemake.output['transformed_df'], sep='\t')

loadings_df = pd.DataFrame(pca.components_.T*np.sqrt(pca.explained_variance_), columns=['x', 'y'], index=df.index)
loadings_df = loadings_df.join(lookup_table)
loadings_df.to_csv(snakemake.output['loadings_df'], sep='\t')

def plot_loadings(loadings_df):
    fig = go.Figure()
    omega_loadings_df = loadings_df[['x','y','stress_type', 'stress_level']].groupby(['stress_type', 'stress_level']).mean()
    omega_loadings_df.index = [f"omega {' '.join(i)}" for i in omega_loadings_df.index.to_flat_index()]
    fig.add_trace(
        go.Scatter(
            x=omega_loadings_df['x'].values,
            y=omega_loadings_df['y'].values,
            name="Omegas",
            legendgroup="Omegas",
            legendgrouptitle_text="Omegas",
            marker_color='black',
            marker_symbol='circle',
            mode="text+markers",
            text=[l.replace('omega', '&#937;') for l in omega_loadings_df.index.values],
            textposition="top center"
        )
    )

    color_dict = {t:c for t, c in zip(loadings_df.stress_type.unique(), px.colors.qualitative.Plotly)}
    marker_dict = {t:c for t, c in zip(loadings_df.stress_level.unique(), ['circle', 'square', 'diamond', 'cross', 'x'])}
    for cycog, row in loadings_df.iterrows():
        group=f"{row['stress_level']} {row['stress_type']}"
        fig.add_trace(
            go.Scatter(
                x=[0, row['x']],
                y=[0, row['y']],
                name=row['gene_name'],
                legendgroup=group,
                legendgrouptitle_text=group,
                marker_color=color_dict[row['stress_type']],
                marker_symbol=marker_dict[row['stress_level']],
                text=f"cycog: {cycog}<br>"+"<br>".join([f"{k}: {v}" for k,v in row.to_dict().items()])
            )
        )

    fig.update_layout(
        xaxis_title=f"PCA1 ({pca.explained_variance_ratio_[0]:.0%} of variance explained)",
        yaxis_title=f"PCA2 ({pca.explained_variance_ratio_[1]:.0%} of variance explained)",
    )
    return fig

fig = plot_loadings(loadings_df)
fig.write_html(snakemake.output['loadings_fig'])