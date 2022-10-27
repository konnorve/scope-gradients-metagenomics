import pandas as pd
import plotly.express as px

transformed_df = pd.read_table(snakemake.input[0])
px.scatter(
    transformed_df,
    x='x',
    y='y',
    color=snakemake.wildcards['column']
).write_image(snakemake.output[0])