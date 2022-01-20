from dataclasses import dataclass
import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import obone


@dataclass
class ALZanalysis:
    directory: str

    def __post_init__(self):
        os.chdir(self.directory)

        self.gene_weights_1 = self._get_json_weights("alz_gene_weights_1.json")

    def _get_json_weights(self, json_file: str):
        with open(json_file) as file_in:
            gene_weights = json.load(file_in)
        gene_weights = {int(k): v for k, v in gene_weights.items()}
        return gene_weights

    def avrampou2019(self) -> obone.BoNE:
        gse138024 = obone.GEO(accessionID="GSE138024")
        survival = gse138024.survival()
        expr = pd.read_parquet("GSE138024-GPL17021-expr.parquet.gzip")
        expr = expr.set_index("Name")
        avrampou = obone.BoNE(expr, survival)

        name = "c genotype"
        groups = ["WT", "RGS4 KO"]
        avrampou.init(name, self.gene_weights_1, groups)
        return avrampou

    def rodriguez2021(self) -> obone.BoNE:
        gse164788 = obone.GEO(accessionID="GSE164788")
        survival = gse164788.survival()
        name = "c drug_concentration_um"
        survival[name] = survival[name].fillna("Control")
        survival[name] = survival[name].astype(str)
        # rename 1 -> 1.0 so that samples 10 and 1 have different regex matches
        survival[name] = survival[name].replace("1", "1.0")
        expr = pd.read_parquet("GSE164788-GPL18573-expr.parquet.gzip")
        expr = expr.set_index("gene_name")
        rodriguez = obone.BoNE(expr, survival)

        name = "c drug_concentration_um"
        groups = ["Control", "1.0", "10"]  # "0.3", "3"]
        rodriguez.init(name, self.gene_weights_1, groups)
        return rodriguez

    def dong2013(self):
        gse40060 = obone.GEO(accessionID="GSE40060")
        survival = gse40060.survival()
        expr = gse40060.expr(rename_genes=True)
        dong = obone.BoNE(expr, survival)

        name = "c source_name_ch1"
        groups = ["endogenous", "overexpressed"]
        dong.init(name, self.gene_weights_1, groups)
        return dong

    def gse1691687(self):
        # raw_files
        pass

    def tan2020(self):
        gse150696 = obone.GEO(accessionID="GSE150696")
        survival = gse150696.survival()
        expr = gse150696.expr()

    def ryan2021(self):
        gse169687 = obone.GEO(accessionID="GSE169687")
        survival = gse169687.survival()
        expr = gse169687.expr()


if __name__ == "__main__":
    import sys

    dir = sys.argv[1]
    alz = ALZanalysis(dir)

    # avrampou = alz.avrampou2019()
    # plt.figure(figsize=(10, 5), dpi=100)
    # avrampou.violin()
    # plt.savefig("avrampou2019_violin.png")

    rodriguez = alz.rodriguez2021()
    # plt.figure(figsize=(10, 5), dpi=100)
    # rodriguez.violin()
    # plt.savefig("rodriguez2020_violin.png")

    # dong = alz.dong2013()
    # plt.figure(figsize=(10, 5), dpi=100)
    # dong.violin()
    # plt.savefig("dong2013_violin.png")
