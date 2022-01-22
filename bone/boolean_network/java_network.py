import subprocess
import os
from typing import Any
import pandas as pd

# java version: jdk1.8.0_45


def stepminer(
    bv,
    file_rl: str = "network.rl",
    p_thr: float = 0.1,
    s_thr: float = 3.0,
    d_thr: float = 0.05,
) -> str:
    if isinstance(bv, pd.DataFrame):
        bv["BitVector"] = bv.apply(lambda x: "".join(x.astype(str)), axis=1)
        bv_file = "bitvector.txt"
        bv["BitVector"].to_csv(bv_file, sep="\t")
    elif isinstance(bv, str):
        bv_file = bv

    file_path = os.path.dirname(os.path.abspath(__file__))
    stepminer = os.path.join(file_path, "stepminer-1.1.jar")
    subprocess.run(
        [
            "java",
            "-cp",
            stepminer,
            "-Xms64m",
            "-Xmx10G",
            "tools.CustomAnalysis",
            "boolean",
            "bitMatrix",
            file_rl,
            bv_file,
            "false_filler.ph",
            "All",
            p_thr,
            s_thr,
            d_thr,
        ]
    )
    os.remove(bv_file)
    subprocess.run(
        [
            "java",
            "-cp",
            stepminer,
            "-Xms64m",
            "-Xmx10G",
            "tools.CustomAnalysis",
            "boolean",
            "bitMatrixFill",
            file_rl,
        ]
    )
    subprocess.run(
        [
            "java",
            "-cp",
            stepminer,
            "-Xms64m",
            "-Xmx10G",
            "tools.CustomAnalysis",
            "boolean",
            "bitMatrixFillStats",
            file_rl,
        ]
    )
    # create -res.txt for human readable results
    file_res_txt = file_rl.split(".")[0] + "-res.txt"
    subprocess.run(
        [
            "java",
            "-cp",
            stepminer,
            "-Xms64m",
            "-Xmx10G",
            "tools.CustomAnalysis",
            "boolean",
            "bitMatrixPrint",
            file_rl,
            ">",
            file_res_txt,
        ]
    )
    return file_rl
