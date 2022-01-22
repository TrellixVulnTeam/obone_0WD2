import pandas as pd
import numpy as np
import scanpy as sc
import tarfile
import os
import re


def normalize(expr: pd.DataFrame, probe_type: str = "cpm") -> pd.DataFrame:
    if probe_type.lower() != "cpm":
        raise ValueError(f"{probe_type} normalization not available. Please use 'cpm'")

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


def add_probeID(expr: pd.DataFrame, organism: str, probe_type: str):
    file_path = os.path.dirname(os.path.abspath(__file__))
    homo_sapiens = os.path.join(file_path, "references/homo_sapiens.csv")
    mus_musculus = os.path.join(file_path, "references/mus_musculus.csv")

    organism = organism.title()
    if organism not in ["Homo Sapiens", "Mus Musculus"]:
        raise ValueError("Only 'Homo Sapiens' and 'Mus Musculus' probeIDs available")

    probe_type = probe_type.upper()
    if probe_type not in ["ENST", "ENSG", "ENSMUST", "ENSMUSG"]:
        raise ValueError("Only 'ENST', 'ENSG', 'ENSMUST', or 'ENSMUSG' are available.")

    if organism == "Homo Sapiens":
        probe_df = pd.read_csv(homo_sapiens, index_col=0)[probe_type]
        expr = expr.merge(probe_df, how="left", right_index=True, left_index=True)
    if organism == "Mus Musculus":
        probe_df = pd.read_csv(mus_musculus, index_col=0)[probe_type]
        expr = expr.merge(probe_df, how="left", right_index=True, left_index=True)

    expr = expr.rename(columns={probe_type: "ProbeID"})
    expr = expr.reset_index()
    expr = expr.set_index(["ProbeID", "Name"])
    return expr
