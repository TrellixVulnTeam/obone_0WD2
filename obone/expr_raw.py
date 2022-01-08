import glob
import tarfile
import os
import shutil
import re
import numpy as np
import pandas as pd


sep = "\t"
skiprows = 1


os.chdir("ALZ")

tar_file = "GSE169687_RAW.tar"
tar_dir = tar_file.split(".")[0]
if not os.path.exists(tar_dir):
    with tarfile.open(tar_file) as tar:
        tar.extractall(tar_dir)

for file in os.listdir(tar_dir):
    gsm = re.sub(".*(GSM[0-9]+).*", "\\1", file)
    filepath = os.path.join(tar_dir, file)
    gsm_df = pd.read_csv(filepath, sep="\t", skiprows=1, usecols=[0, 6])
    try:
        gsm_df.columns = ["ProbeID", gsm]
    except:
        raise ValueError("gsm_df should only contain ID column and read data")
    exclude = [
        "no_feature",
        "ambiguous",
        "too_low_aQual",
        "not_aligned",
        "alignment_not_unique",
    ]
    gsm_df = gsm_df[~gsm_df["ProbeID"].isin(exclude)]
    gsm_df = gsm_df.set_index("ProbeID")
    if "df" not in locals():
        df = gsm_df
    else:
        df = df.merge(gsm_df, left_index=True, right_index=True)

shutil.rmtree(tar_dir)
# # Add gene names
# scaff_path = "/booleanfs2/sahoo/Data/SeqData/genome/"
# if gsm_df["ProbeID"].str.contains("ENSG").any():
#     gene_names = pd.read_csv(
#         scaff_path + "Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.t.txt",
#         sep="\t",
#         names=["T", "ProbeID", "Name"],
#     )
#     del gene_names["T"]
#     gsm_df["ProbeID"] = gsm_df["ProbeID"].str.split(".", expand=True)[0]
#     gsm_df = gsm_df.merge(gene_names, how="left", on="ProbeID")

# elif gsm_df["ProbeID"].str.contains("ENST").any():
#     gene_names = pd.read_csv(
#         scaff_path + "Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.t.txt",
#         sep="\t",
#         names=["ProbeID", "G", "Name"],
#     )
#     del gene_names["G"]
#     gsm_df["ProbeID"] = gsm_df["ProbeID"].str.split(".", expand=True)[0]
#     gsm_df = gsm_df.merge(gene_names, how="left", on="ProbeID")

# elif gsm_df["ProbeID"].str.contains("ENSMUST").any():
#     gene_names = pd.read_csv(
#         scaff_path + "Mus_musculus.GRCm38.94.chr_patch_hapl_scaff.t.txt",
#         sep="\t",
#         names=["ProbeID", "G", "Name"],
#     )
#     del gene_names["G"]
#     gsm_df["ProbeID"] = gsm_df["ProbeID"].str.split(".", expand=True)[0]
#     gsm_df = gsm_df.merge(gene_names, how="left", on="ProbeID")

# elif gsm_df["ProbeID"].str.contains("ENSMUSG").any():
#     gene_names = pd.read_csv(
#         scaff_path + "Mus_musculus.GRCm38.94.chr_patch_hapl_scaff.t.txt",
#         sep="\t",
#         names=["T", "ProbeID", "Name"],
#     )
#     del gene_names["T"]
#     gsm_df["ProbeID"] = gsm_df["ProbeID"].str.split(".", expand=True)[0]
#     gsm_df = gsm_df.merge(gene_names, how="left", on="ProbeID")

# else:
#     gsm_df["Name"] = gsm_df["ProbeID"]

# # gsm_df = gsm_df.merge(gene_names, how='left', on='ProbeID') Move Name column to correct spot
# names = gsm_df["Name"]
# del gsm_df["Name"]
# gsm_df.insert(1, "Name", names)

# gsm_df = gsm_df.replace(np.NaN, 0)
