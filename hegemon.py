from dataclasses import dataclass
import pandas as pd
import numpy as np
from io import BytesIO
import warnings


@dataclass
class Hegemon:
    expr: pd.DataFrame
    survival: pd.DataFrame

    def __post_init__(self):
        # check for na in expr
        if self.expr.isnull().values.any():
            raise ValueError("Expr contains null values. Please remove.")

        # change expr index names
        if self.expr.index.name != "Name" and "Name" not in self.expr.index.names:
            raise ValueError("Expr index must contain column 'Name'")

        # change survival index name to Sample and make all capital letters
        self.survival.index.name = "Sample"
        self.survival.index = self.survival.index.str.upper()

    def thr(self) -> pd.DataFrame:
        warnings.filterwarnings("ignore")
        array = np.array(self.expr)
        array = np.sort(array)

        mean = np.mean(array, axis=1).reshape(-1, 1)
        sstot = np.subtract(array, mean)
        sstot = np.sum(np.square(sstot), axis=1)
        count1 = np.array(range(1, array.shape[1] + 1))
        count2 = np.array(range(array.shape[1] - 1, -1, -1))

        sum1 = np.cumsum(array, axis=1)
        sum1 = np.insert(sum1, 0, 0, axis=1)
        sum1 = sum1[:, :-1]
        sum2 = np.sum(array, axis=1)
        sum2 = sum2.reshape(-1, 1) - array.cumsum(axis=1)
        val = np.sum(array, axis=1).reshape(-1, 1)
        sum2 = np.concatenate((val, sum2), axis=1)
        sum2 = sum2[:, :-1]

        m1 = sum1 / (count1 - 1)
        m1[:, 0] = 0
        m2 = sum2 / (count2 + 1)
        m2[:, 0] = np.squeeze(mean)

        tmp1 = array + sum1
        tmp1 = tmp1 / count1
        tmp1 = mean - tmp1
        tmp1 = np.square(tmp1)
        tmp1 = tmp1 * count1
        sum1sq = array - mean
        sum1sq = np.square(sum1sq)
        sum1sq = sum1sq - tmp1
        to_add = (count1 - 1) * np.square(mean - m1)
        sum1sq = sum1sq + to_add
        sum1sq = np.cumsum(sum1sq, axis=1)

        tmp2 = sum2 - array
        tmp2 = tmp2 / count2
        tmp2 = mean - tmp2
        tmp2 = np.square(tmp2)
        tmp2 = tmp2 * count2
        to_sub = np.negative(np.square(array - mean))
        to_add = (count2 + 1) * np.square(mean - m2)
        to_sub = to_sub + to_add - tmp2
        to_sub = np.cumsum(to_sub, axis=1)
        sum2sq = np.dstack([sstot] * array.shape[1])
        sum2sq = np.squeeze(sum2sq)
        sum2sq = sum2sq + to_sub

        sse = sum1sq + sum2sq
        sse[:, -1] = sstot

        min_index = np.argmin(sse, axis=1)
        m1 = np.array([np.mean(row[: i + 1]) for row, i in zip(array, min_index)])
        m2 = np.array([np.mean(row[i + 1 :]) for row, i in zip(array, min_index)])

        thr = (m1 + m2) / 2
        thr = thr.reshape(-1, 1)

        min_sse = np.min(sse, axis=1)
        min_sse[min_sse < 1e-10] = 0
        statistic = sstot - min_sse
        if array.shape[1] <= 4:
            statistic = statistic / min_sse / 2
        elif array.shape[1] > 4:
            to_div = min_sse / (array.shape[1] - 4)
            statistic = statistic / to_div
            statistic = statistic / 3
        statistic = np.where(min_sse == 0, 0, statistic)
        statistic = statistic.reshape(-1, 1)

        thr = np.concatenate((thr, statistic), axis=1)
        thr = pd.DataFrame(
            thr, columns=["Threshold", "Statistic"], index=self.expr.index
        )
        thr["Low Noise"] = thr["Threshold"] - 0.5
        thr["High Noise"] = thr["Threshold"] + 0.5
        warnings.filterwarnings("default")
        return thr

    def bv(self) -> pd.DataFrame:
        bv = np.array(self.expr)
        thr = self.thr()
        low = np.array(thr["Low Noise"]).reshape(-1, 1)
        high = np.array(thr["High Noise"]).reshape(-1, 1)
        bv = np.where(bv < low, 0, bv)
        bv = np.where(bv > high, 2, bv)
        between = np.logical_and(bv >= low, bv <= high)
        bv = np.where(between, 1, bv)
        bv = pd.DataFrame(bv, columns=self.expr.columns, index=self.expr.index)
        bv = bv.astype(int)
        return bv

    def idx(self) -> pd.DataFrame:
        idx = self.expr.copy()
        df_bin = BytesIO()
        idx.to_csv(df_bin, mode="wb", sep="\t")
        df_bin.seek(0)

        ptr = []
        pos = 0
        for line in df_bin:
            if pos == 0:
                pos += len(line)
            else:
                ptr.append(pos)
                pos += len(line)

        idx["Ptr"] = ptr
        return idx["Ptr"]

    def ih(self, name: str) -> pd.DataFrame:
        ih_df = self.survival[name]
        ih_df.insert(1, "ArrayHeader", ih_df["Sample"])
        ih_df.insert(1, "ClinicalHeader", ih_df["c title"])
        return ih_df
