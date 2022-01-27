import GEOparse
import pandas as pd
import numpy as np
import re
import os
import glob
import warnings

from dataclasses import dataclass

from .preprocess import *


@dataclass
class GEO:
    accessionID: str
    remove_soft: bool = False

    def init(self):
        # Add GEOparse gsms and gpls attributes as class attributes
        gse = GEOparse.get_GEO(geo=self.accessionID, silent=True)
        self.gsms = gse.gsms
        self.gpls = gse.gpls
        self.default_gpl = list(self.gpls.keys())[0]

        # remove downloaded soft file
        if self.remove_soft:
            os.remove(glob.glob(f"{self.accessionID}*.soft*")[0])

    def survival(self, gpl_name: str = None) -> pd.DataFrame:
        """Creates metadata information for each GSM (sample)

        Args:
            gpl_name (string): gpl name associated

        Returns:
            pd.DataFrame: survival dataframe including all samples and metadata
        """
        if not hasattr(self, "gsms"):
            self.init()

        # use first gpl if none specified
        if gpl_name == None:
            gpl_name = self.default_gpl
            print(f"Survival: Using default GPL: {gpl_name}")

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
        norm_type: str = None,
        log2: bool = False,
        get_genes: bool = True,
        get_gsm: str = None,
    ) -> pd.DataFrame:
        """Create expression file from class attribute (NCBI GEO Accession ID)

        Args:
            gpl_name (str, optional): NCBI GEO GPL. Defaults to None.
            norm_type (str, optional): Normalization method. Defaults to None.
            log2 (bool, optional): Take a log2 of entire expression file. Defaults to False.
            get_genes (bool, optional): rename genes from GEOparse gpl table. Defaults to True.
            get_gsm (str, optional): Rename sample (columns) with survival column. Defaults to None.

        Raises:
            ValueError: if GEOparse gsm.table is empty

        Returns:
            pd.DataFrame: Expression file.
        """
        if not hasattr(self, "gsms"):
            self.init()

        if gpl_name == None:
            gpl_name = self.default_gpl
            print(f"Expr: Using default GPL: {gpl_name}")

        for name, gsm in self.gsms.items():
            # confirm that gsm correlates to called gpl
            if gsm.metadata["platform_id"][0] != gpl_name:
                continue

            # check if table contains expr data
            if gsm.table.empty:
                raise ValueError("soft file does not contain expr data")

            # collect expr data
            gsm_df = gsm.table[["ID_REF", "VALUE"]]
            gsm_df.columns = ["ProbeID", name]
            gsm_s = gsm_df.set_index("ProbeID")

            # merge each column into one dataframe
            if "expr" not in locals():
                expr = gsm_s
            else:
                expr = expr.merge(gsm_s, left_index=True, right_index=True)

        if norm_type:
            expr = normalize(expr, norm_type)

        if log2:
            expr = np.log2(expr + 1)

        if get_genes:
            expr = self._get_genes(expr, gpl_name)

        if get_gsm:
            expr = self._get_gsm(expr, gpl_name, get_gsm)

        if expr.isnull().values.any():
            warnings.warn(
                "Expression file contains NA values. Must be filled prior to use with BoNE"
            )
        return expr

    def _get_genes(self, expr: pd.DataFrame, gpl_name: str):
        # clean original gpl table column names
        gpl_table = self.gpls[gpl_name].table
        rename = lambda x: x.replace("_", " ").title().replace(" ", "")
        gpl_table.columns = [rename(c) for c in gpl_table.columns]
        gpl_table = gpl_table.rename(columns={"Id": "ProbeID"})

        # rename columns containing symbol info to 'Symbol'
        sym_rename = []
        [sym_rename.append(c) for c in gpl_table.columns if "Sym" in c]
        sym_rename = {r: "Name" for r in sym_rename}
        gpl_table = gpl_table.rename(columns=sym_rename)
        gpl_table["Name"] = gpl_table["Name"].replace("---", pd.NA)
        gpl_table["Name"] = gpl_table["Name"].str.upper()

        # rename columns containing description info to 'Description'
        def_rename = ["GeneName", "Definition"]
        [def_rename.append(c) for c in gpl_table.columns if "Desc" in c]
        def_rename = {r: "Description" for r in def_rename}
        gpl_table = gpl_table.rename(columns=def_rename)

        # merge with expr frame to rename with gene names
        gpl_table = gpl_table[["ProbeID", "Name", "Description"]]
        expr = expr.merge(gpl_table, how="left", on="ProbeID")
        expr = expr.set_index(["ProbeID", "Name"])

        # drop description for now, perhaps include in future version
        expr = expr.drop("Description", axis=1)
        return expr

    def _get_gsm(self, expr: pd.DataFrame, gpl_name: str, survival_col: str):
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
