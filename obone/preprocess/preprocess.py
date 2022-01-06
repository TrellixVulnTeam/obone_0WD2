import pandas as pd
import scanpy as sc
import tarfile
import glob
import gzip
import pathlib
import shutil
import os


def counts2expr(counts_df: pd.DataFrame) -> pd.DataFrame:
    counts_df = counts_df.set_index(["ProbeID", "Name"])
    adata = sc.AnnData(counts_df.T)
    # adata.var_names_make_unique()
    # 1e6 = counts per million (cpm) normalization
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata, base=2)
    return adata.to_df().T


def tar2rcc(tar_file: str) -> str:
    """Extracts all files from a tarfile into a new directory with '_RCC' appended

    Args:
        tar_file (str): must be a tarfile

    Returns:
        list: Path to new '_RCC' directory
    """
    if not tarfile.is_tarfile(tar_file):
        raise ValueError(f"{tar_file} is not a tarfile")

    filename = os.path.basename(os.path.abspath(tar_file))
    filename = filename.split(".")[0]
    rcc_dir = filename + "_RCC"

    with tarfile.open(tar_file) as tar:
        tar.extractall(rcc_dir)

    gunzip_dir(rcc_dir)
    return rcc_dir


def gunzip_dir(my_dir: str) -> None:
    """Unzips any gzipped files in given directory

    Args:
        my_dir (str): path to directory for files to be gunzipped
    """
    for gzip_file in glob.glob(os.path.join(my_dir, "*.gz")):
        gunzip_file = pathlib.Path(gzip_file).stem
        gunzip_file = os.path.join(my_dir, gunzip_file)
        with gzip.open(gzip_file, "rb") as file_in, open(gunzip_file, "wb") as file_out:
            shutil.copyfileobj(file_in, file_out)
            os.remove(gzip_file)
