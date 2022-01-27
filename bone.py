from dataclasses import dataclass
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib.colors as colors
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns

from .hegemon import Hegemon


@dataclass
class BoNE(Hegemon):
    def __post_init__(self):
        super().__post_init__()

    def init(self, survival_col: str, gene_weights: dict, groups: dict) -> None:
        self.survival_col = survival_col
        self.gene_weights = self._clean_gene_weights(gene_weights)
        self.groups = self._clean_groups(groups)
        self.sample_data = self.data()

    def _clean_gene_weights(self, gene_weights) -> dict:
        # check self.gene_weights is dict
        if not isinstance(gene_weights, dict):
            raise ValueError(
                "'gene_weights' must be provided as dict: {\weight: [genes]}"
            )
        # guarantee that gene weight is float
        gene_weights = {float(weight): genes for weight, genes in gene_weights.items()}
        return gene_weights

    def _clean_groups(self, groups) -> dict:
        # make group a dictionary
        if not isinstance(groups, dict):
            if not isinstance(groups, list):
                groups = list(groups)
            groups = {value: [value] for value in groups}

        # ensure group values (samples) are a list
        for key in groups.keys():
            value = groups[key]
            if not isinstance(value, list):
                groups[key] = [value]

        # check that no group values encapsulate other values for later regex
        group_keys = list(groups.keys())
        for i, val in enumerate(group_keys):
            w_o_val = group_keys[:i] + group_keys[i + 1 :]
            for v in w_o_val:
                if val in v:
                    raise ValueError(
                        "Group values cannot contain other values (i.e. in [1, 10], 10 contains 1. Change '1' to '1.0')"
                    )
        return groups

    def _check_init(self) -> None:
        if not hasattr(self, "survival_col"):
            raise AttributeError(
                "BoNE must first be initialized for this method. Please use the '.init()' method"
            )

    def score(self) -> pd.DataFrame:
        self._check_init()

        for weight, group in self.gene_weights.items():
            expr = self.expr
            expr = expr[expr.index.get_level_values("Name").isin(group)]
            thr = self.thr()
            thr = thr[thr.index.get_level_values("Name").isin(group)]

            sd = np.std(expr, axis=1).replace(0, 1)
            v = expr.sub(thr.iloc[:, 3], axis=0)
            v = v.div(3)
            v = v.div(sd, axis=0)

            rank = v.sum(axis=0)
            rank.name = f"Weight: {weight}"

            if "ranks" not in locals():
                ranks = rank
            else:
                ranks = pd.concat([ranks, rank], axis=1)

        weights = [float(w) for w in list(self.gene_weights.keys())]
        ranks["Score"] = np.dot(ranks, np.array(weights))
        ranks.name = "Sample"
        ranks = ranks["Score"]

        survival = self.survival[self.survival_col]
        survival.name = "Sample Type"

        score_df = pd.concat((survival, ranks), join="inner", axis=1)
        score_df.index.name = "Sample"
        return score_df

    def data(self) -> pd.DataFrame:
        df = self.score()

        # map group_val to samples and groups
        all_sample_types = []
        group_val_group = {}
        group_val_sample_type = {}
        for i, (group_name, sample_types) in enumerate(self.groups.items()):
            # craete single list of all samples in all groups
            all_sample_types.extend(sample_types)
            # map group_val to group name for later use
            group_val_group[i] = group_name
            # create group_val per sample type
            for sample_type in sample_types:
                group_val_sample_type[sample_type] = i
        # set-up regex to filter dataframe
        searchfor = "|".join(all_sample_types)
        df = df[df["Sample Type"].str.contains(searchfor)]
        df["Sample Type"] = df["Sample Type"].str.extract(f"({searchfor})")
        df["group_val"] = df["Sample Type"].replace(group_val_sample_type)

        # add Group annotation
        df = df.reset_index()
        df["Group"] = df.groupby(["group_val"])["Sample"].transform("count")
        df["Group"] = "(" + df["Group"].astype(str) + ")"
        df["Group"] = df["group_val"].replace(group_val_group).str.title() + df["Group"]

        # add color
        color = {i: get_cmap("Paired")(i) for i in df["group_val"].unique()}
        df["group_color"] = df["group_val"].map(color)

        # add pvalue and roc_auc score
        df["Pval"] = None
        df["ROC AUC"] = None
        control_scores = list(df[df["group_val"] == 0]["Score"])
        for val in range(1, df["group_val"].max() + 1):
            # pval
            group_score = list(df[df["group_val"] == val]["Score"])
            _, pval = ttest_ind(control_scores, group_score, equal_var=False)
            if pval < 0.05:
                df.loc[df["group_val"] == val, "Pval"] = pval

            # roc score
            roc_auc_data = df[df["group_val"].isin([0, val])]
            roc_auc = roc_auc_score(roc_auc_data["group_val"], roc_auc_data["Score"])
            df.loc[df["group_val"] == val, "ROC AUC"] = roc_auc

        # sort data by group_val for proper coloring
        df = df.sort_values("group_val")
        return df

    def title_bar(self):
        ax = plt.subplot2grid((4, 1), (0, 0))
        group_vals = list(self.sample_data.sort_values("Score")["group_val"])
        axis = [0, len(group_vals), 0, 5]
        ax.axis(axis)
        cmap = colors.ListedColormap(self.sample_data["group_color"])

        ax.imshow(
            np.array(group_vals).reshape(1, -1),
            interpolation="nearest",
            cmap=cmap,
            extent=axis,
            aspect="auto",
        )
        ax.set(xticks=range(len(group_vals)), xticklabels=[], yticklabels=[])
        ax.tick_params(top=False, left=False, bottom=False, right=False)
        ax.grid(which="major", color="black", alpha=0.5, linestyle="-", linewidth=0.75)

        res = set(self.sample_data["ROC AUC"].dropna())
        if len(res) > 0:
            res_text = f'ROC: {",".join([str(round(val,2)) for val in res])}'
            ax.text(len(group_vals), 4, res_text)

        # make title_bar top patches and annotation
        divider = make_axes_locatable(ax)
        ax1 = divider.append_axes("top", size="100%", pad="20%", frame_on=False)
        ax1.axison = False
        ax1.axis([0, len(group_vals), 0, 5])
        ax1.grid(False)

        spacer = len(group_vals) / (max(group_vals) + 1)
        for i, index in enumerate(self.sample_data.drop_duplicates("group_val").index):
            ax1.add_patch(
                patches.Rectangle(
                    (i * spacer, 0),
                    1,
                    3,
                    facecolor=self.sample_data["group_color"][index],
                    alpha=1.0,
                )
            )
            ax1.text(
                i * spacer + 1,
                1,
                self.sample_data["Group"][index],
                rotation="horizontal",
                ha="left",
                va="center",
                fontsize=12,
            )
        return ax

    def violin(self):
        self._check_init()

        ax = self.title_bar()
        ax = plt.subplot2grid((4, 1), (1, 0), rowspan=3)
        sns.set_theme(palette=set(self.sample_data["group_color"]))

        ax = sns.violinplot(
            x="Score",
            y="Group",
            data=self.sample_data,
            inner="quartile",
            linewidth=0.5,
            ax=ax,
        )
        ax = sns.swarmplot(
            x="Score",
            y="Group",
            color="blue",
            alpha=0.2,
            ax=ax,
            data=self.sample_data,
        )
        ax.set(ylabel=None)
        ax.xaxis.grid(True, clip_on=False)

        if "Pval" in self.sample_data.columns:
            text = self.sample_data.drop_duplicates("group_val")
            text = text[text["Pval"].notnull()]
            y_value = 0.5
            for group in text.index:
                ax.text(
                    text.loc[group, "Score"],
                    y_value,
                    f'p={text.loc[group, "Pval"]:.1e}',
                    horizontalalignment="center",
                    size=12,
                    color="0.3",
                )
                y_value += 1
        return ax

    def density(self):
        self._check_init()

        ax = self.title_bar()
        ax = plt.subplot2grid((4, 1), (1, 0), rowspan=3)

        for i in self.sample_data.drop_duplicates("group_val").index:
            group = self.sample_data.loc[i, "Group"]
            df = self.sample_data[self.sample_data["Group"] == group]
            s = df.reset_index()["index"]
            if len(s) != 1:
                ax = s.plot.kde(
                    bw_method=1.0,
                    ax=ax,
                    c=self.sample_data["group_color"][i],
                    label=group,
                )
            elif len(s) == 1:
                df["y"] = 1
                ax = df.plot.line(
                    x=group,
                    y="y",
                    ax=ax,
                    c=self.sample_data["group_color"][i],
                    label=group,
                )
                ax.axvline(x=df.index[0], c=self.sample_data["group_color"][i])
