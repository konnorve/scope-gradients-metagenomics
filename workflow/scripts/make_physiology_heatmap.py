import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
from pathlib import Path

all_phys_data = pd.read_table(snakemake.input[0], index_col=['cruise', 'sample_id'])
out_dir = Path(snakemake.output[0])
out_dir.mkdir()

for cruise in all_phys_data.index.get_level_values('cruise').unique():
    phys_data = all_phys_data.loc[cruise].sort_values(
            by=snakemake.config['physiology heatmap']['sort order'],
            ascending=snakemake.config['physiology heatmap']['ascending']
    )
    phys_data = phys_data[snakemake.config['physiology heatmap']['colors'].keys()].astype(float)
    norm_phys_data=(phys_data-phys_data.min())/(phys_data.max() - phys_data.min())

    img = np.zeros((len(phys_data), len(snakemake.config['physiology heatmap']['colors']), 4))

    for j, (col, cmap) in enumerate(snakemake.config['physiology heatmap']['colors'].items()):
        norm = norm_phys_data[col].to_numpy()
        cmap = cm.get_cmap(cmap)
        img[:,j,:] = cmap(norm)

    plt.figure(figsize = (20,5))
    plt.imshow(img, aspect='auto')
    plt.yticks([])
    plt.xticks(range(len(snakemake.config['physiology heatmap']['colors'])), snakemake.config['physiology heatmap']['colors'].keys(), rotation = 45)
    plt.tight_layout()
    plt.savefig(out_dir / f"{cruise}_physiology_heatmap.png")
