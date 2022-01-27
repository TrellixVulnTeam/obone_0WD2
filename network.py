import pandas as pd
from io import StringIO
import re
import tempfile
import subprocess
import os
import sys
from typing import Any

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
            # create temporary bitvector file
            with tempfile.NamedTemporaryFile(delete=False) as temp:
                self.bv["BitVector"].to_csv(temp.name, sep="\t")
                self.bv = temp.name

        file_path = os.path.dirname(os.path.abspath(__file__))
        self.stepminer_path = os.path.join(file_path, "references/stepminer-1.1.jar")

        self._step1()
        self._step2()
        self._step3()
        self.network = self._step4()
        os.remove(self.file_rl)

    def _step1(self):
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
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
        )

    def _step2(self):
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
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
        )

    def _step3(self):
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
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
        )

    def _step4(self):
        output = subprocess.run(
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
            ],
            capture_output=True,
            text=True,
        )
        output = output.stdout
        pattern = "AID\tnorel\tlohi\tlolo\thihi\thilo\teqv\topp"
        span_start = re.search(pattern, output).span()[0]
        output = output[span_start:]
        df = pd.read_csv(StringIO(output), sep="\t", index_col=0)
        df.index.name = "ProbeID"
        return df
