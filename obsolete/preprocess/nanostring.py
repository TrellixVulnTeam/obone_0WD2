from dataclasses import dataclass
import pandas as pd
import os
from bs4 import BeautifulSoup
from io import StringIO


@dataclass
class NanoString:
    rcc_files: list

    def __post_init__(self):
        # map sample name to parsed rcc file with {attribute: pd.DataFrame}
        sample_hash = {}
        for rcc_file in self.rcc_files:
            sample = self.parse_rcc(rcc_file)
            sample_name = os.path.basename(rcc_file)
            sample_hash[sample_name] = sample

            # make all sample attribute class attribute if does not exist
            for attribute in sample.keys():
                if not hasattr(self, attribute):
                    setattr(self, attribute, None)

        # compile all RCC attribute into class attribute dataframe
        for attribute in list(vars(self)):
            if attribute != "rcc_files":
                for sample_name, sample in sample_hash.items():
                    to_merge = sample[attribute]
                    to_merge.columns = [sample_name]
                    if "df" not in locals():
                        df = to_merge
                    else:
                        df = df.merge(to_merge, left_index=True, right_index=True)
                setattr(self, attribute, df)
                del df

    def parse_rcc(self, file_in):
        with open(file_in) as file_in:
            soup = BeautifulSoup(file_in.read(), "html.parser")
        tags = [tag.name for tag in soup.find_all()]

        sample = {}
        for tag in tags:
            tag_string = getattr(soup, tag).string.strip()
            if tag == "code_summary":
                df = pd.read_csv(StringIO(tag_string))
                df = df[df["Name"].notnull()]
                df = df.set_index(["CodeClass", "Name", "Accession"])
                df = df.astype("float")
            else:
                df = pd.read_csv(StringIO(tag_string), names=["Attribute", "Value"])
                df = df.set_index("Attribute")
            df.columns.name = tag
            sample[tag] = df
        return sample
