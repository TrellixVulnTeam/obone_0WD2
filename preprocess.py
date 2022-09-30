import pandas as pd
import numpy as np
import scanpy as sc
import tempfile
from pathlib import Path
import tarfile
import os
import re


def log2(expr: pd.DataFrame):
    expr = np.log2(expr + 1)
    return expr


def log_normalize(expr: pd.DataFrame, norm_type: str = "cpm") -> pd.DataFrame:
    """Normalize given expression file via scanpy preprocessing

    Args:
        expr (pd.DataFrame): Expression file to be normalized
        norm_type (str, optional): Normalization method. Defaults to "cpm".
        log2 (bool, optional): Take log2 of all values after normalization. Defaults to False.

    Raises:
        ValueError: Only cpm normalization is currently available

    Returns:
        pd.DataFrame: Normalized expression dataframe
    """
    if norm_type.lower() != "cpm":
        raise ValueError(f"{norm_type} normalization not available. Please use 'cpm'")
    elif norm_type.lower == "cpm":
        # 1e6 = counts per million (cpm) normalization
        norm_type = 1e6

    adata = sc.AnnData(expr.T)
    sc.pp.normalize_total(adata, target_sum=norm_type, exclude_highly_expressed=True)
    sc.pp.log1p(adata, base=2)
    expr = adata.to_df().T
    return expr


def read_raw(tar_file: str, **kwargs) -> pd.DataFrame:
    """Parse tar file into expression frame

    Args:
        tar_file (str): path to tarfile

    Raises:
        ValueError: **kwargs must clean tarfile to a two column frame with ID and expression value

    Returns:
        pd.DataFrame: Expression dataframe
    """    
    if "sep" not in kwargs:
        kwargs["sep"] = "\t"

    tar_dir = tempfile.TemporaryDirectory()
    with tarfile.open(tar_file) as tar:
        def is_within_directory(directory, target):
            
            abs_directory = os.path.abspath(directory)
            abs_target = os.path.abspath(target)
        
            prefix = os.path.commonprefix([abs_directory, abs_target])
            
            return prefix == abs_directory
        
        def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
        
            for member in tar.getmembers():
                member_path = os.path.join(path, member.name)
                if not is_within_directory(path, member_path):
                    raise Exception("Attempted Path Traversal in Tar File")
        
            tar.extractall(path, members, numeric_owner) 
            
        
        safe_extract(tar, tar_dir.name)
    for file in os.listdir(tar_dir.name):
        gsm = re.sub(".*(GSM[0-9]+).*", "\\1", file)
        file = os.path.join(tar_dir.name, file)
        gsm_df = pd.read_csv(file, **kwargs)
        try:
            gsm_df.columns = ["ID", gsm]
        except:
            print(gsm_df)
            raise ValueError(
                "gsm_df should only contain ID column and read data. Enter **kwargs for pd.read_csv() clean"
            )
        gsm_df = gsm_df.set_index("ID")
        if "df" not in locals():
            df = gsm_df
        else:
            df = df.merge(gsm_df, left_index=True, right_index=True)
    return df


def add_probeID(expr: pd.DataFrame, probe_type: str) -> pd.DataFrame:
    """Only RNAseq. Add an index level of ProbeID. ProbeID will be merged with existing index

    Args:
        expr (pd.DataFrame): Expression file to parse
        organism (str): Type of organism for reference probe
        probe_type (str): Probe type to use for ProbeID values

    Raises:
        ValueError: Only certain probe types are currently available

    Returns:
        pd.DataFrame: Frame with ProbeID in index
    """
    file_path = os.path.dirname(os.path.abspath(__file__))
    homo_sapiens = os.path.join(file_path, "references/homo_sapiens.csv")
    mus_musculus = os.path.join(file_path, "references/mus_musculus.csv")

    probe_type = probe_type.upper()
    if probe_type not in ["ENST", "ENSG", "ENSMUST", "ENSMUSG"]:
        raise ValueError("Only 'ENST', 'ENSG', 'ENSMUST', or 'ENSMUSG' are available.")
    elif probe_type in ["ENST", "ENSG"]:
        probe_df = pd.read_csv(homo_sapiens, index_col=0)[probe_type]
        probe_df = probe_df[~probe_df.index.duplicated(keep="first")]
        expr = expr.merge(probe_df, how="left", right_index=True, left_index=True)
    elif probe_type in ["ENSMUST", "ENSMUSG"]:
        probe_df = pd.read_csv(mus_musculus, index_col=0)[probe_type]
        probe_df = probe_df[~probe_df.index.duplicated(keep="first")]
        expr = expr.merge(probe_df, how="left", right_index=True, left_index=True)

    expr = expr.rename(columns={probe_type: "ProbeID"})
    if expr.index.name:
        current_index = [expr.index.name]
    elif expr.index.names:
        current_index = expr.index.names
    expr = expr.reset_index()
    expr = expr.set_index(["ProbeID"] + current_index)
    return expr
