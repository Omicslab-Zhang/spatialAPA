

import os
import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm

# ---------- 解析参数 ----------
parser = argparse.ArgumentParser(description="APA coordinate remapping")
parser.add_argument("--sample_info", required=True, help="Path to sample_info.txt")
parser.add_argument("--window", type=int, default=25, help="Window size for merging APA coordinates")
parser.add_argument("--outdir", default=".", help="Directory to save outputs")
args = parser.parse_args()

'''
usage:
python /public/home/jiangzh/projects/sc_apa/APA_dataset_combine/spatialAPA_merge_pas.py \
    --sample_info /public/home/jiangzh/projects/sc_apa/APA_dataset_combine/test_apamerge.txt \
    --window 25 \
    --outdir /public/home/jiangzh/projects/sc_apa/APA_dataset_combine/test_tmp_merge

'''

# ---------- 加载 sample_info ----------
sample_info = pd.read_csv(args.sample_info, sep="\t", header=None, names=["sample", "filepath"])
window = args.window
outdir = args.outdir

# ---------- 统一 APA 坐标格式 ----------
def normalize_coord(coord):
    if "-" in coord.split(":")[1]:  # chr1:1311566-1311689:-
        chrom, pos_range, strand = coord.split(":")
        end = pos_range.split("-")[1]
        return f"{chrom}:{end}:{strand}"
    parts = coord.split(":")
    ### scape input
    if len(parts) == 4:  # chr1:1311566:27:-
        return f"{parts[0]}:{parts[1]}:{parts[3]}"
    return coord  # already standard: chr1:1311566:-

# ---------- 解析并收集所有 APA 坐标 ----------
all_coords = []
sample_matrices = {}

# def read_expression_file(fpath):
#     ext = os.path.splitext(fpath)[1].lower()
#     if ext == '.csv':
#         return pd.read_csv(fpath, index_col=0)
#     elif ext in ['.txt', '.tsv']:
#         return pd.read_csv(fpath, sep='\t', index_col=0)
#     else:
#         raise ValueError(f"Unsupported file format: {fpath}")

def read_expression_file(fpath):
    fpath = str(fpath)
    if fpath.endswith('.csv') or fpath.endswith('.csv.gz'):
        return pd.read_csv(fpath, index_col=0)
    elif fpath.endswith('.tsv') or fpath.endswith('.tsv.gz') or fpath.endswith('.txt') or fpath.endswith('.txt.gz'):
        return pd.read_csv(fpath, sep='\t', index_col=0)
    else:
        raise ValueError(f"Unsupported file format: {fpath}")


for _, row in tqdm(sample_info.iterrows(), total=len(sample_info), desc="Reading samples"):
    sample_id = row["sample"]
    fpath = row["filepath"]
    if not os.path.exists(fpath):
        print(f"Warning: File not found: {fpath}")
        continue
    #df = pd.read_csv(fpath, index_col=0)
    df = read_expression_file(fpath)
    df.index = df.index.map(normalize_coord)
    sample_matrices[sample_id] = df
    all_coords.extend(df.index.tolist())

# ---------- 构建坐标 DataFrame ----------
def parse_coords(coords):
    parts = [x.split(":") for x in coords]
    return pd.DataFrame({
        'chr': [p[0] for p in parts],
        'pos': [int(p[1]) for p in parts],
        'strand': [p[2] for p in parts],
        'old_coord': coords
    })

coord_df = parse_coords(list(set(all_coords)))

# ---------- 坐标合并 ----------
def merge_coords_with_mapping(df, window=25):
    df = df.sort_values(['chr', 'strand', 'pos']).reset_index(drop=True)
    mapping = []
    i = 0
    while i < len(df):
        group = [df.loc[i]]
        j = i + 1
        while j < len(df):
            row = df.loc[j]
            if row['chr'] != group[0]['chr'] or row['strand'] != group[0]['strand']:
                break
            if all(abs(row['pos'] - g['pos']) <= window for g in group):
                group.append(row)
                j += 1
            else:
                break
        merged_pos = min(r['pos'] for r in group)
        merged_coord = f"{group[0]['chr']}:{merged_pos}:{group[0]['strand']}"
        for r in group:
            mapping.append({
                'old_coord': r['old_coord'],
                'new_coord': merged_coord
            })
        i = j
    return pd.DataFrame(mapping)

mapping_df = merge_coords_with_mapping(coord_df, window=window)

# ---------- 重映射每个样本 ----------
os.makedirs(outdir, exist_ok=True)

for sample_id, mat in tqdm(sample_matrices.items(), desc="Remapping"):
    mapping = mapping_df[mapping_df['old_coord'].isin(mat.index)].drop_duplicates('old_coord')
    mat = mat.copy()
    mat['new_coord'] = mat.index.map(mapping.set_index('old_coord')['new_coord'])
    mat = mat.dropna(subset=['new_coord'])
    remapped = mat.groupby('new_coord').sum()
    # 输出
    sample_dir = os.path.join(outdir, sample_id)
    os.makedirs(sample_dir, exist_ok=True)
    remapped.to_csv(os.path.join(sample_dir, "remapped_apa.csv"))
    mapping.to_csv(os.path.join(sample_dir, "apa_mapping.csv"), index=False)
    print(f"{sample_id}: {mat.shape[0]} → {remapped.shape[0]} (rows)")
