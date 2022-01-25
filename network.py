import subprocess
import os
from typing import Any
import pandas as pd

from dataclasses import dataclass

from .bone import BoNE
from .geo import GEO

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

        self.step1()
        self.step2()
        self.step3()
        self.human_readable()

    def setp1(self):
        """
        export CLASSPATH="./stepminer-1.1.jar"
        stepminer="java -cp $CLASSPATH -Xms64m -Xmx10G tools.CustomAnalysis"
        ${stepminer} boolean bitMatrix $FILE.rl \
            data/peters-2017-ibd-bv.txt \
            $FILE.ph All 0.1 2.5 0.05
        """
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

    def step2(self):
        """
        export CLASSPATH="./stepminer-1.1.jar"
        stepminer="java -cp $CLASSPATH -Xms64m -Xmx10G tools.CustomAnalysis"
        cp $FILE.rl $FILE-1.rl
        ${stepminer} boolean bitMatrixFill $FILE-1.rl
        """
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

    def step3(self):
        """
        export CLASSPATH="./stepminer-1.1.jar"
        stepminer="java -cp $CLASSPATH -Xms64m -Xmx10G tools.CustomAnalysis"
        ${stepminer} boolean bitMatrixFillStats $FILE-1.rl
        """
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

    def human_readable(self):
        """
        export CLASSPATH="./stepminer-1.1.jar"
        stepminer="java -cp $CLASSPATH -Xms64m -Xmx10G tools.CustomAnalysis"
        ${stepminer} boolean bitMatrixPrint $FILE-1.rl > $FILE-res.txt
        """
        # create -res.txt for human readable results
        file_res_txt = self.file_rl.split(".")[0] + "-res.txt"
        subprocess.run(
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
            ]
        )


if __name__ == "__main__":
    gse40060 = GEO(accessionID="GSE40060")
    survival = gse40060.survival()
    expr = gse40060.expr(rename_genes=True, probeID="ENSG")
    expr = expr.fillna(0)
    dong = BoNE(expr, survival)
    bv = dong.bv()
    Stepminer(bv)
