import pandas as pd
import logging
from pathlib import Path

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

def read_type(s):
    if "_1_unmerged" in s:
        return "unmerged"
    if "_2_unmerged" in s:
        return "unmerged_2"
    if "_merged" in s:
        return "merged"
    if "_unpaired" in s:
        return "unpaired"

read_ser_list = []
base_ser_list = []
for k, v in snakemake.input.items():
    logging.info(f"{k}\t{v}")
    df = pd.read_table(v)

    df['id'] = df.file.apply(lambda s: str(Path(s).name).split('.')[0])
    df['sample'] = df.id.apply(lambda s: s.replace("_1_unmerged", "").replace("_2_unmerged", "").replace("_merged", "").replace("_unpaired", ""))
    df['read_type'] = df.id.apply(lambda s: read_type(s))

    assert df[df.read_type == "unmerged"].num_seqs.to_list() == df[df.read_type == "unmerged_2"].num_seqs.to_list()
    
    base_ser = df.groupby(['sample']).sum().sum_len
    base_ser.name = k
    base_ser_list.append(base_ser)

    df = df[df.read_type != "unmerged_2"]
    read_ser = df.groupby(['sample']).sum().num_seqs
    read_ser.name = k
    read_ser_list.append(read_ser)

read_df = pd.concat(read_ser_list, axis=1)
read_df = read_df.sort_index()
read_df = read_df.transpose()
read_df.to_csv(snakemake.output['read_count'], sep='\t')

base_df = pd.concat(base_ser_list, axis=1)
base_df = base_df.sort_index()
base_df = base_df.transpose()
base_df.to_csv(snakemake.output['base_count'], sep='\t')