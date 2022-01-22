from io import StringIO
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
            self.bv = bv_file

        file_path = os.path.dirname(os.path.abspath(__file__))
        self.stepminer_path = os.path.join(file_path, "references/stepminer-1.1.jar")

    def build(self):
        subprocess.run(
            [
                "java",
                "-cp",
                self.stepminer_path,
                "-Xms64m",
                "-Xmx10G",
                "tools.CustomAnalysis",
                "boolean",
                "bitMatrix",
                self.file_rl,
                self.bv,
                "false_filler.ph",
                "All",
                str(self.p_thr),
                str(self.s_thr),
                str(self.d_thr),
            ]
        )
        subprocess.run(
            [
                "java",
                "-cp",
                self.stepminer_path,
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
                self.stepminer_path,
                "-Xms64m",
                "-Xmx10G",
                "tools.CustomAnalysis",
                "boolean",
                "bitMatrixFillStats",
                self.file_rl,
            ]
        )

    def get_readable(self):
        # create -res.txt for human readable results
        file_res_txt = self.file_rl.split(".")[0] + "-res.txt"
        network_res = subprocess.run(
            [
                "java",
                "-cp",
                self.stepminer_path,
                "-Xms64m",
                "-Xmx10G",
                "tools.CustomAnalysis",
                "boolean",
                "bitMatrixPrint",
                self.file_rl,
                ">",
                file_res_txt,
            ],
            capture_output=True,
        )
        df = pd.read_csv(StringIO(network_res), sep="\t")
        df.to_csv(file_res_txt)
        return df
