import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import sys

from dataclasses import dataclass

sys.path.insert(0, "/Users/oliver.tucher/Code/obone")
import bone


@dataclass
class ALZanalysis:
    def __post_init__(self):
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        self.gene_weights_1 = self._get_json_weights("gene_weights_1.json")

    def _get_json_weights(self, json_file: str):
        with open(json_file) as file_in:
            gene_weights = json.load(file_in)
        gene_weights = {int(k): v for k, v in gene_weights.items()}
        return gene_weights

    def avrampou2019(self) -> bone.BoNE:
        gse138024 = bone.GEO(accessionID="GSE138024")
        survival = gse138024.survival()
        expr = pd.read_parquet("GSE138024-GPL17021-expr.parquet.gzip")
        expr = expr.set_index("Name")
        avrampou = bone.BoNE(expr, survival)

        name = "c genotype"
        groups = ["WT", "RGS4 KO"]
        avrampou.init(name, self.gene_weights_1, groups)
        return avrampou

    def rodriguez2021_xyz(self) -> bone.BoNE:
        gse164788 = bone.GEO(accessionID="GSE164788")

        survival = gse164788.survival()
        name = "c drug_concentration_um"
        survival[name] = survival[name].fillna("Control")
        survival[name] = survival[name].astype(str)
        # rename 1 -> 1.0 so that samples 10 and 1 have different regex matches
        survival[name] = survival[name].replace("1", "1.0")
        # filter survival for specific group
        survival = survival[survival["c title"].str.contains("dsRNA + lipofectamine")]

        expr = pd.read_parquet("GSE164788-GPL18573-expr.parquet.gzip")
        expr = expr.set_index("gene_name")
        rodriguez = bone.BoNE(expr, survival)

        name = "c drug_concentration_um"
        groups = ["Control", "1.0", "10"]  # "0.3", "3"]
        rodriguez.init(name, self.gene_weights_1, groups)
        return rodriguez

    def dong2013(self):
        gse40060 = bone.GEO(accessionID="GSE40060")
        survival = gse40060.survival()
        expr = gse40060.expr(rename_genes=True)
        expr = expr.fillna(0)
        expr = bone.add_probeID(expr, "Homo Sapiens", "ENSG")
        dong = bone.BoNE(expr, survival)

        name = "c source_name_ch1"
        groups = ["endogenous", "overexpressed"]
        dong.init(name, self.gene_weights_1, groups)
        return dong

    def gse1691687(self):
        # raw_files
        pass

    def tan2020(self):
        gse150696 = bone.GEO(accessionID="GSE150696")
        survival = gse150696.survival()
        expr = gse150696.expr()

    def ryan2021(self):
        gse169687 = bone.GEO(accessionID="GSE169687")
        survival = gse169687.survival()
        print(survival)

    def peters2017(self):
        gse83687 = bone.GEO(accessionID="GSE83687")
        survival = gse83687.survival()

        expr_file = f"GSE83687-{gse83687.default_gpl}-expr.parquet.gzip"
        if not os.path.exists(expr_file):
            expr = bone.read_raw("GSE83687_RAW.tar")
            expr.index = expr.index.str.upper()
            expr.index.name = "Name"
            expr = bone.add_probeID(expr, "Homo Sapiens", "ENSG")
            expr = bone.normalize(expr, log2=True)
            expr.to_parquet(expr_file, compression="gzip")
        else:
            expr = pd.read_parquet(expr_file)

        peters = bone.BoNE(expr, survival)
        bv = peters.bv()
        bone.Stepminer(bv)


def violin(bone_obj):
    plt.figure(figsize=(10, 5), dpi=100)
    bone_obj.violin()
    plt.savefig("violin.png")


if __name__ == "__main__":
    alz = ALZanalysis()
    peters = alz.peters2017()
