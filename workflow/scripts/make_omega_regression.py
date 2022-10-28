import pandas as pd
import plotly.express as px
from pathlib import Path

nutrient_concentrations = pd.read_table(snakemake.input['physiology_data'], index_col="sample_id")

omegas = pd.read_table(snakemake.input['omegas'])
out_dir = Path(snakemake.output[0])
out_dir.mkdir()

omegas = omegas.melt(id_vars='sample_id', var_name='stress', value_name='omega')
df = omegas.join(nutrient_concentrations, on='sample_id')
df['N_P_Ratio'] = df['NO3_NO2'] / df['PO4']
df = df.melt(id_vars=['sample_id', 'stress', 'omega'], value_vars=['N_P_Ratio','SiO4','NO3_NO2','TN','TON','PO4','TP','TOP','TOC','NH4','NO2'], var_name='nut_type', value_name='nut_conc')

px.scatter(
    df[(df['stress'].str.contains('Phosphorus')) & (df['nut_type'].isin(['PO4','TP','TOP']))], 
    x='nut_conc', 
    y='omega', 
    facet_col='nut_type', 
    color='stress',
    width=1800, height=600,
    log_x=True, trendline="ols", trendline_options={'log_x':True}
).write_image(out_dir / f"phosphorus.png")

for nut_type in ['PO4','TP','TOP']:
    px.scatter(
        df[(df['stress'].str.contains('Phosphorus')) & (df['nut_type']==nut_type)], 
        x='nut_conc', 
        y='omega', 
        facet_col='nut_type', 
        color='stress',
        log_x=True, trendline="ols", trendline_options={'log_x':True}
    ).write_image(out_dir / f"{nut_type}_phosphorus_stress.png")

px.scatter(
    df[(df['stress'].str.contains('Nitrogen')) & (df['nut_type'].isin(['N_P_Ratio','NO3_NO2','TN','TON','NH4','NO2']))], 
    x='nut_conc', 
    y='omega', 
    facet_col='nut_type', 
    color='stress',
    width=3600, height=600,
    log_x=True, trendline="ols", trendline_options={'log_x':True}
).write_image(out_dir / f"nitrogen.png")

for nut_type in ['N_P_Ratio','NO3_NO2','TN','TON','NH4','NO2']:
    px.scatter(
        df[(df['stress'].str.contains('Nitrogen')) & (df['nut_type']==nut_type)], 
        x='nut_conc', 
        y='omega', 
        facet_col='nut_type', 
        color='stress',
        log_x=True, trendline="ols", trendline_options={'log_x':True}
    ).write_image(out_dir / f"{nut_type}_nitrogen_stress.png")