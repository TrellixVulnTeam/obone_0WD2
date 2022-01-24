from dataclasses import dataclass
import GEOparse
import pandas as pd
import numpy as np
import re
import scanpy as sc
import os
import glob


@dataclass
class GEO:
    accessionID: str

    def _geo_init(self):
        # Add GEOparse gsms and gpls attributes as class attributes
        gse = GEOparse.get_GEO(geo=self.accessionID, silent=True)
        self.gsms = gse.gsms
        self.gpls = gse.gpls
        self.default_gpl = list(self.gpls.keys())[0]

        # remove downloaded soft file
        os.remove(glob.glob(f"{self.accessionID}*.soft*")[0])

    def survival(self, gpl_name: str = None) -> pd.DataFrame:
        """Creates metadata information for each GSM (sample)

        Args:
            gpl_name (string): gpl name associated

        Returns:
            pd.DataFrame: survival dataframe including all samples and metadata
        """
        if not hasattr(self, "gsms"):
            self._geo_init()

        # use first gpl if none specified
        if gpl_name == None:
            gpl_name = self.default_gpl
            print("Survival: Using default GPL")

        to_drop = [
            "geo_accession",
            "status",
            "date",
            "protocol",
            "proccesing",
            "data_processing",
            "contact",
            "supplementary",
            "platform_id",
            "series_id",
            "relation",
            "channel",
        ]

        all_metadata = {}
        for name, gsm in self.gsms.items():
            # confirm gsm is associated with defined gpl
            gsm_gpl = gsm.metadata["platform_id"][0]
            if gsm_gpl != gpl_name:
                continue

            # remove keys from metadata that aren't desired in survival
            metadata = gsm.metadata.copy()
            for key in gsm.metadata:
                for drop in to_drop:
                    if re.search(drop, key.lower()):
                        metadata.pop(key, None)
            all_metadata[name] = metadata
        all_metadata = pd.DataFrame(all_metadata).T

        # split columns with multiple values into seperate columns
        # convert all list values into string to split
        all_metadata = all_metadata.applymap(lambda x: "\t".join(x), na_action="ignore")
        for column in all_metadata.columns:
            to_merge = all_metadata[column].str.split("\t", expand=True)
            if len(to_merge.columns) > 1:
                # create column names for additional columns
                col_names = [str(i + 1) for i in range(len(to_merge.columns))]
                col_names = [column + "_" + name for name in col_names]
                to_merge.columns = col_names
            else:
                to_merge.columns = [column]

            if "df" not in locals():
                df = to_merge
            else:
                df = df.merge(to_merge, left_index=True, right_index=True)

        df.index.name = "Sample"

        # rename columns with value
        for column in df.columns:
            if df[column].str.contains(":").all():
                value = df[column].str.extract(r"(.*:)").iloc[0, 0]
                value = str(value).strip(":")
                df = df.rename(columns={column: value})
                if isinstance(df[value], pd.DataFrame):
                    continue
                # remove 'xxxx: ' from each column value
                df[value] = df[value].str.extract(r".*: (.*)")

        # add 'c ' to all columns per hegemon requirements
        df = df.rename(columns={col: "c " + col for col in df.columns})
        return df

    def expr(
        self,
        gpl_name: str = None,
        log2: bool = False,
        normalize: bool = False,
        rename_genes: bool = False,
        rename_samples_with: str = None,
    ) -> pd.DataFrame:
        if not hasattr(self, "gsms"):
            self._geo_init()

        if gpl_name == None:
            gpl_name = self.default_gpl
            print("Expr: Using default GPL")

        for name, gsm in self.gsms.items():
            # confirm that gsm correlates to called gpl
            if gsm.metadata["platform_id"][0] != gpl_name:
                continue

            # check if table contains expr data
            if gsm.table.empty:
                raise ValueError("soft file does not contain expr data")

            # collect expr data
            gsm_df = gsm.table
            gsm_df.columns = ["Name", name]
            gsm_s = gsm_df.set_index("Name")

            if log2:
                gsm_s = np.log2(gsm_s + 1)

            # merge each column into one dataframe
            if "expr" not in locals():
                expr = gsm_s
            else:
                expr = expr.merge(gsm_s, left_index=True, right_index=True)

        if normalize:
            expr = self.normalize(expr, "cpm")

        if rename_genes:
            expr = self.rename_genes(expr, gpl_name)

        if rename_samples_with != None:
            expr = self.rename_samples_with(expr, gpl_name, rename_samples_with)

        return expr

    def rename_genes(self, expr: pd.DataFrame, gpl_name: str):
        # clean original gpl table column names
        gpl_table = self.gpls[gpl_name].table
        rename = lambda x: x.replace("_", " ").title().replace(" ", "")
        gpl_table.columns = [rename(c) for c in gpl_table.columns]
        gpl_table = gpl_table.rename(columns={"Id": "Name"})

        # rename columns containing symbol info to 'Symbol'
        sym_rename = []
        [sym_rename.append(c) for c in gpl_table.columns if "Sym" in c]
        sym_rename = {r: "Symbol" for r in sym_rename}
        gpl_table = gpl_table.rename(columns=sym_rename)
        gpl_table["Symbol"] = gpl_table["Symbol"].replace("---", pd.NA)
        gpl_table["Symbol"] = gpl_table["Symbol"].str.upper()

        # rename columns containing description info to 'Description'
        def_rename = ["GeneName", "Definition"]
        [def_rename.append(c) for c in gpl_table.columns if "Desc" in c]
        def_rename = {r: "Description" for r in def_rename}
        gpl_table = gpl_table.rename(columns=def_rename)

        # merge with expr frame to rename with gene names
        gpl_table = gpl_table[["Name", "Symbol", "Description"]]
        expr = expr.merge(gpl_table, how="left", on="Name")
        expr.loc[expr["Symbol"].isna(), "Symbol"] = expr["Name"]
        expr = expr.drop("Name", axis=1).rename(columns={"Symbol": "Name"})
        expr = expr.set_index("Name")

        # drop description for now, perhaps include in future version
        expr = expr.drop("Description", axis=1)
        return expr

    def rename_samples_with(self, expr: pd.DataFrame, gpl_name: str, survival_col: str):
        survival = self.survival(gpl_name)

        if survival_col not in survival.columns:
            raise ValueError(f"{survival_col} not in Survival frame")

        survival = survival[survival_col]
        survival = survival.to_dict()
        rename_dict = {v: k for k, v in survival.items()}
        expr = expr.rename(columns=rename_dict)

        columns = list(rename_dict.keys())
        not_in_expr = [c for c in columns if c not in expr.columns]
        if len(not_in_expr) > 0:
            print(f"{not_in_expr} columns were not renamed.")
        return expr

    def normalize(expr: pd.DataFrame, type: str = "cpm") -> pd.DataFrame:
        if type.lower() != "cpm":
            raise ValueError(f"{type} normalization not available. Please use 'cpm'")

        adata = sc.AnnData(expr.T)
        # target_sum of 1e6 equates to counts per million (cpm) normalization
        if type == "cpm":
            target_sum = 1e6
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata, base=2)
        expr = adata.to_df().T
        return expr
