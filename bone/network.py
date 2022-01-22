import subprocess
import os
from typing import Any
import pandas as pd

from dataclasses import dataclass

# java version: jdk1.8.0_45


@dataclass
class Stepminer:
    bv: Any
    file_rl: str = "network.rl"
    p_thr: float = 0.1
    s_thr: float = 3.0
    d_thr: float = 0.05

    def __post_init__(self):
        if isinstance(self.bv, pd.DataFrame):
            self.bv["BitVector"] = self.bv.apply(
                lambda x: "".join(x.astype(str)), axis=1
            )
            bv_file = "bitvector.txt"
            self.bv["BitVector"].to_csv(bv_file, sep="\t")
        elif isinstance(self.bv, str):
            bv_file = self.bv

        file_path = os.path.dirname(os.path.abspath(__file__))
        stepminer = os.path.join(file_path, "references/stepminer-1.1.jar")

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
                self.file_rl,
                bv_file,
                "false_filler.ph",
                "All",
                str(self.p_thr),
                str(self.s_thr),
                str(self.d_thr),
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
                self.file_rl,
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
                self.file_rl,
            ]
        )
        # create -res.txt for human readable results
        file_res_txt = self.file_rl.split(".")[0] + "-res.txt"
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
                self.file_rl,
                ">",
                file_res_txt,
            ]
        )
