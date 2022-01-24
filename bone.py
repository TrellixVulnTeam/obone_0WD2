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
        self.sample_data = self.plot_data(survival_col, gene_weights, groups)
        self.cval = list(self.sample_data.sort_values("Score")["Cval"])
        self.cval_colors = list(self.sample_data["Color"])

        group_data = self.sample_data.drop_duplicates("Annotation")
        group_data = group_data.set_index("Annotation")
        self.group_data = group_data
        self.annotations = list(group_data.index)
        self.group_colors = list(group_data["Color"])
        self.res = list(group_data["ROC AUC"].dropna())

    def rank(self, gene_weights: dict) -> pd.DataFrame:
        for weight, group in gene_weights.items():
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

        weights = [float(w) for w in list(gene_weights.keys())]
        ranks["Score"] = np.dot(ranks, np.array(weights))
        ranks.name = "Sample"
        return ranks

    def score(self, survival_col: str, gene_weights: dict) -> pd.DataFrame:
        survival = self.survival[survival_col]
        survival.name = "Sample Type"

        # add score
        score = self.rank(gene_weights)["Score"]
        # FUTURE: raise warning if survival index != expr columns and alert how many rows lost
        df = pd.concat((survival, score), join="inner", axis=1)
        df.index.name = "Sample"
        return df

    def plot_data(self, survival_col, gene_weights, groups) -> pd.DataFrame:
        df = self.score(survival_col, gene_weights)

        # make group a dictionary
        if not isinstance(groups, dict):
            groups = {value: [value] for value in groups}
        group_keys = list(groups.keys())

        # check that no group values encapsulate other values for later regex
        for i, val in enumerate(group_keys):
            w_o_val = group_keys[:i] + group_keys[i + 1 :]
            for v in w_o_val:
                if val in v:
                    raise ValueError(
                        "Group values cannot contain other values (i.e. in [1, 10], 10 contains 1. Change '1' to '1.0')"
                    )

        # map cval to samples and groups
        all_sample_types = []
        cval_group = {}
        cval_sample_type = {}
        for i, (group_name, sample_types) in enumerate(groups.items()):
            # ensure samples are a list
            if not isinstance(sample_types, list):
                sample_types = [sample_types]
            # craete single list of all samples in all groups
            all_sample_types.extend(sample_types)
            # map cval to group name for later use
            cval_group[i] = group_name
            # create cval per sample type
            for sample_type in sample_types:
                cval_sample_type[sample_type] = i
        searchfor = "|".join(all_sample_types)
        df = df[df["Sample Type"].str.contains(searchfor)]
        df["Sample Type"] = df["Sample Type"].str.extract(f"({searchfor})")
        df["Cval"] = df["Sample Type"].replace(cval_sample_type)

        # add annotation
        df = df.reset_index()
        df["Annotation"] = df.groupby(["Cval"])["Sample"].transform("count")
        df["Annotation"] = "(" + df["Annotation"].astype(str) + ")"
        df["Annotation"] = df["Cval"].replace(cval_group).str.title() + df["Annotation"]
        df = df.set_index("Sample")

        # add color
        color = {i: get_cmap("Paired")(i) for i in range(len(groups.keys()))}
        df["Color"] = df["Cval"].map(color)

        # add pvalue and roc_auc score
        control_scores = list(df[df["Cval"] == 0]["Score"])
        for val in df[df["Cval"] != 0]["Cval"].unique():
            group_score = list(df[df["Cval"] == val]["Score"])
            _, pval = ttest_ind(control_scores, group_score, equal_var=False)
            if pval < 0.05:
                df.loc[df["Cval"] == val, "Pval"] = pval

            roc_auc_data = df[df["Cval"].isin([0, val])]
            roc_auc = roc_auc_score(roc_auc_data["Cval"], roc_auc_data["Score"])
            df.loc[df["Cval"] == val, "ROC AUC"] = roc_auc

        # sort data by cval for proper coloring
        df = df.sort_values("Cval")
        return df

    def title_bar(self):
        ax = plt.subplot2grid((4, 1), (0, 0))
        cval = np.array(self.cval).reshape(1, len(self.cval))
        extent = [0, len(self.cval), 0, 5]
        ax.axis(extent)
        cmap = colors.ListedColormap(self.cval_colors)

        ax.imshow(
            cval,
            interpolation="nearest",
            cmap=cmap,
            extent=extent,
            aspect="auto",
        )
        ax.set(xticks=range(len(self.cval)), xticklabels=[], yticklabels=[])
        ax.tick_params(top=False, left=False, bottom=False, right=False)
        for _, spine in ax.spines.items():
            spine.set_visible(False)
        ax.grid(which="major", color="black", alpha=0.5, linestyle="-", linewidth=0.75)

        res_text = f'ROC: {",".join([str(round(val,2)) for val in self.res])}'
        ax.text(len(self.cval), 4, res_text)
        return ax

    def title_bar_top(self, ax):
        divider = make_axes_locatable(ax)
        ax1 = divider.append_axes("top", size="100%", pad="20%", frame_on=False)
        ax1.axison = False
        ax1.axis([0, len(self.cval), 0, 5])
        ax1.grid(False)

        spacer = len(self.cval) / len(self.annotations)
        for i in range(len(self.annotations)):
            ax1.add_patch(
                patches.Rectangle(
                    (i * spacer, 0),
                    1,
                    3,
                    facecolor=self.group_colors[i],
                    edgecolor="none",
                    alpha=1.0,
                )
            )
            ax1.text(
                i * spacer + 1,
                1,
                self.annotations[i],
                rotation="horizontal",
                ha="left",
                va="center",
                fontsize=12,
            )

    def violin(self, survival_col=None, gene_weights=None, groups=None):
        if not hasattr(self, "sample_data"):
            if survival_col == None:
                raise ValueError("Either init() BoNE or pass arguments")
            self.init(survival_col, gene_weights, groups)

        ax = self.title_bar()
        self.title_bar_top(ax)

        ax = plt.subplot2grid((4, 1), (1, 0), rowspan=3)
        sns.set_theme(palette=self.group_colors)

        ax = sns.violinplot(
            x="Score",
            y="Annotation",
            data=self.sample_data,
            inner="quartile",
            linewidth=0.5,
            ax=ax,
        )
        ax = sns.swarmplot(
            x="Score",
            y="Annotation",
            color="blue",
            alpha=0.2,
            ax=ax,
            data=self.sample_data,
        )
        ax.set(ylabel=None)
        ax.xaxis.grid(True, clip_on=False)

        if "Pval" in self.sample_data.columns:
            text = self.group_data[self.group_data["Pval"].notnull()]
            y_value = 0.5
            for annotation in text.index:
                ax.text(
                    text.loc[annotation, "Score"],
                    y_value,
                    f'p={text.loc[annotation, "Pval"]:.1e}',
                    horizontalalignment="center",
                    size=12,
                    color="0.3",
                )
                y_value += 1
        return ax

    def density(self, survival_col=None, gene_weights=None, groups=None):
        if not hasattr(self, "sample_data"):
            if survival_col == None:
                raise ValueError("Either init() BoNE or pass arguments")
            self.init(survival_col, gene_weights, groups)

        ax = self.title_bar()
        self.title_bar_top(ax)

        ax = plt.subplot2grid((4, 1), (1, 0), rowspan=3)
        df = self.sample_data.reset_index(drop=True)
        for i in range(len(self.annotations)):
            df1 = df[df["Annotation"] == self.annotations[i]]
            annotation = df1["Annotation"].iloc[0]
            s = df1.reset_index()["index"]
            s.name = annotation
            if len(s) != 1:
                ax = s.plot.kde(
                    bw_method=1.0, ax=ax, c=self.group_colors[i], label=annotation
                )
            elif len(s) == 1:
                df1["y"] = 1
                ax = df1.plot.line(
                    x=annotation, y="y", ax=ax, c=self.group_colors[i], label=annotation
                )
                ax.axvline(x=df1.index[0], c=self.group_colors[i])
