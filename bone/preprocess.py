import pandas as pd
import numpy as np
import scanpy as sc


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
