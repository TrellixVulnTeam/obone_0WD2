import pandas as pd
import numpy as np
import scanpy as sc
import tarfile
import os
import re


def normalize(expr: pd.DataFrame, type: str = "cpm") -> pd.DataFrame:
    if type.lower() != "cpm":
        raise ValueError(f"{type} normalization not available. Please use 'cpm'")

    expr = expr.set_index(["ProbeID", "Name"])
    adata = sc.AnnData(expr.T)
    # adata.var_names_make_unique()
    # 1e6 = counts per million (cpm) normalization
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata, base=2)
    expr = adata.to_df().T
    return expr


def log2(expr: pd.DataFrame):
    expr = np.log2(expr + 1)
    return expr


def read_raw(tar_file, **kwargs):
    tar_dir = tar_file.split(".")[0]
    if not os.path.exists(tar_dir):
        with tarfile.open(tar_file) as tar:
            tar.extractall(tar_dir)

    for file in os.listdir(tar_dir):
        gsm = re.sub(".*(GSM[0-9]+).*", "\\1", file)
        filepath = os.path.join(tar_dir, file)
        gsm_df = pd.read_csv(filepath, **kwargs)
        try:
            gsm_df.columns = ["ProbeID", gsm]
        except:
            raise ValueError("gsm_df should only contain ID column and read data")
        gsm_df = gsm_df[gsm_df["ProbeID"].str.contains("ENSG|ENST|ENSMUST|ENSMUSG")]
        gsm_df = gsm_df.set_index("ProbeID")
        if "df" not in locals():
            df = gsm_df
        else:
            df = df.merge(gsm_df, left_index=True, right_index=True)
