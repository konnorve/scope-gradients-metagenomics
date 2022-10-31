import pandas as pd
import plotly.express as px
from pathlib import Path

nutrient_concentrations = pd.read_table(snakemake.input['physiology_data'], index_col="sample_id")
omegas = pd.read_table(snakemake.input['omegas'], index_col="sample_id")
out_dir = Path(snakemake.output[0])
out_dir.mkdir()

df = omegas.join(nutrient_concentrations)
df = pd.melt(df, 
    id_vars=nutrient_concentrations.columns, 
    value_vars=omegas.columns, 
    var_name='omega_name', 
    value_name='omega_value')

for nutrient_type, nutrient_concentration in snakemake.config['nutrient regressions'].items():
    nutrient_type_omegas = [omega for omega, info in snakemake.config['omegas'] if info['attributes']['nutrient']==nutrient_type]
    px.scatter(
        df[df['omega_name'].isin(nutrient_type_omegas)], 
        x=nutrient_concentration, 
        y='omega_value', 
        color='omega_name',
        log_x=True, trendline="ols", trendline_options={'log_x':True}
    ).write_image(out_dir / f"{nutrient_type}_{nutrient_concentration}_nitrogen_stress.png")