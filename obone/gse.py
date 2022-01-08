from dataclasses import dataclass
import GEOparse
import pandas as pd
import numpy as np
import re
import os
import glob


@dataclass
class GSE:
    accessionID: str

    def geoparse_init(self):
        # Add GEOparse gsms and gpls attributes as class attributes
        gse = GEOparse.get_GEO(geo=self.accessionID, silent=True)
        self.gsms = gse.gsms
        self.gpls = gse.gpls

        # # remove downloaded soft file
        # os.remove(glob.glob(f"{self.accessionID}*.soft*")[0])

    def survival(self, gpl_name: str = None) -> pd.DataFrame:
        """Creates metadata information for each GSM (sample)

        Args:
            gpl (string): gpl name associated

        Returns:
            pd.DataFrame: survival dataframe including all samples and metadata
        """
        if not hasattr(self, "gsms"):
            self.geoparse_init()
        # use first gpl if none specified
        if gpl_name == None:
            gpl_name = list(self.gpls.keys())[0]
            print(f"Survival GPL: {gpl_name}")

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
                    if re.search(drop, key):
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
                df[value] = df[value].str.extract(r".*: (.*)")

        # add 'c ' label per hegemon requirements
        df = df.rename(columns={col: "c " + col for col in df.columns})

        self._survival = df
        return df

    def expr(self, gpl_name: str, take_log2: bool = False) -> pd.DataFrame:
        if not hasattr(self, "gsms"):
            self.geoparse_init()

        # confirm that gsm correlates to called gpl
        for name, gsm in self.gsms.items():
            if gsm.metadata["platform_id"][0] != gpl_name:
                continue

            # check if table contains expr data
            if gsm.table.empty:
                print("soft file does not contain expr data")
                return

            # collect expr data
            gsm_df = gsm.table.rename(columns={"ID_REF": "Name"})
            gsm_df = gsm_df.set_index("Name")
            gsm_df.columns = [name]

            if take_log2:
                gsm_df[name] = np.where(gsm_df[name] > 0, np.log2(gsm_df[name]), 0)

            # merge each column into one dataframe
            if "expr" not in locals():
                expr = gsm_df
            else:
                expr = expr.merge(gsm_df, left_index=True, right_index=True)

        # expr = self.rename_genes(expr, gpl_name)

        self._expr = expr
        return expr

    def rename_genes(self, expr: pd.DataFrame, gpl_name: str):
        gpl_table = self.gpls[gpl_name].table
        sym_rename = ["GeneSym"]
        [sym_rename.append(c) for c in gpl_table.columns if "Symbol" in c]
        sym_rename = {r: "Symbol" for r in sym_rename}
        gpl_table = gpl_table.rename(columns=sym_rename)

        def_rename = ["GeneName", "GeneDesc", "Definition"]
        [def_rename.append(c) for c in gpl_table.columns if "Description" in c]
        def_rename = {r: "Description" for r in def_rename}
        gpl_table = gpl_table.rename(columns=def_rename)

        gpl_table = gpl_table.rename(columns={"Id": "Name"})
        gpl_table["Symbol"] = gpl_table["Symbol"].replace("---", pd.NA)
        gpl_table = gpl_table[["Name", "Symbol", "Description"]]
        expr = expr.merge(gpl_table, how="left", on="Name")
        expr.loc[expr["Symbol"].isna(), "Symbol"] = expr["Name"]
        expr = expr.drop("Name", axis=1).rename(columns={"Symbol": "Name"})
        expr = expr.set_index("Name")
        expr = expr.drop("Description", axis=1)
        return expr

    def to_gsm(self, expr: pd.DataFrame, survival_col: str):
        if hasattr(self, "_survival"):
            survival = self._survival
        else:
            survival = self.survival()

        survival = survival[survival_col]
        survival = survival.to_dict()
        rename_dict = {v: k for k, v in survival.items()}
        expr = expr.rename(columns=rename_dict)

        columns = list(rename_dict.keys())
        not_in_expr = [c for c in columns if c not in expr.columns]
        if len(not_in_expr) > 0:
            print(f"{not_in_expr} columns were not renamed.")
        return expr
