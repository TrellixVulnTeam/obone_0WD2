import pandas as pd
import numpy as np
import os
import array
import sys
from networkx.utils.union_find import UnionFind


sys.path.insert(0, "/Users/oliver.tucher/Code/obone")


class BINetwork:
    def __init__(self, filename):
        if not os.path.isfile(filename):
            print("Can't open file {0} <br>".format(filename))
            exit()
        self.fp = open(filename, "rb")
        self.file = filename
        self.magic = -1
        self.major = -1
        self.minor = -1
        self.num = -1
        self.numBits = -1
        self.numBytes = -1
        self.startPtr = -1
        self.name = None
        self.low = []
        self.high = []
        self.balanced = []
        self.balancedhash = {}

    def init(self):
        fh = self.fp
        fh.seek(0)
        self.magic = array.array("B", fh.read(1))[0]
        self.major = array.array("B", fh.read(1))[0]
        self.minor = array.array("B", fh.read(1))[0]
        self.low = self.getLow()
        self.high = self.getHigh()
        self.balanced = self.getBalanced()
        self.balancedhash = {}
        for i in range(len(self.balanced)):
            self.balancedhash[self.balanced[i]] = i
        fh.seek(3 + 3 * 4)
        buffer = fh.read(4)
        ptr = array.array("I", buffer)[0]
        buffer = fh.read(4)
        self.num = array.array("I", buffer)[0]
        buffer = fh.read(4)
        self.numBits = array.array("I", buffer)[0]
        self.numBytes = int(self.num * self.numBits / 8) + 1
        fh.seek(ptr)
        self.startPtr = ptr
        return

    def getLow(self):
        return self.readList(0)

    def getHigh(self):
        return self.readList(1)

    def getBalanced(self):
        return self.readList(2)

    def readList(self, num):
        fh = self.fp
        fh.seek(3 + num * 4)
        buffer = fh.read(4)
        ptr = array.array("I", buffer)[0]
        fh.seek(ptr)
        buffer = fh.read(4)
        length = array.array("I", buffer)[0]
        buffer = fh.read(length)
        name = buffer.decode("utf-8")
        buffer = fh.read(4)
        size = array.array("I", buffer)[0]
        res = []
        for i in range(size):
            buffer = fh.read(4)
            length = array.array("I", buffer)[0]
            buffer = fh.read(length)
            res += [buffer.decode("utf-8")]
        return res

    def readBlockID(self, id1):
        if id1 in self.balancedhash:
            a = self.balancedhash[id1]
            return self.readBlockIndex(a)
        return None

    def readBlockIndex(self, a):
        ptr = self.startPtr + 2 * a * self.numBytes
        self.fp.seek(ptr)
        return self.readBlock()

    def readBlock(self):
        fh = self.fp
        buffer1 = fh.read(self.numBytes)
        buffer2 = fh.read(self.numBytes)
        res = []
        for i in range(int(self.num / 2)):
            code = self.readCode(buffer1, buffer2, i, self.numBits)
            res += [code]
        return res

    def readCode(self, buffer1, buffer2, i, numBits, debug=0):
        index = 2 * i * numBits
        v1 = self.readVal(buffer1, index, numBits)
        v2 = self.readVal(buffer1, index + numBits, numBits)
        v3 = self.readVal(buffer2, index, numBits)
        v4 = self.readVal(buffer2, index + numBits, numBits)
        if debug == 1:
            print(index, index + numBits, numBits, v1, v2, v3, v4)
        total = v1 + v2 + v3 + v4
        if total == 1:
            if v1 == 1:
                return 1
            if v2 == 1:
                return 2
            if v3 == 1:
                return 3
            if v4 == 1:
                return 4
        if total == 2:
            if v2 == 1 and v3 == 1:
                return 5
            if v1 == 1 and v4 == 1:
                return 6
        return 0

    def readVal(self, buffer, i, nBits, debug=0):
        index = int(i * nBits / 8)
        v1 = 0
        v2 = 0
        if index < len(buffer):
            v1 = buffer[index]
        if (index + 1) < len(buffer):
            v2 = buffer[index + 1]
        shift = (i * nBits) % 8
        mask = (1 << nBits) - 1
        val = ((v1 | v2 << 8) >> shift) & mask
        if debug == 1:
            print(v1, v2, shift, mask, val)
        return val


class BIGraph:
    def readEqGraph(self, cfile):
        edges = {}
        nodes = {}
        count = 0
        with open(cfile, "r") as netFile:
            for ln in netFile:
                if not ln.startswith("Found : 5"):
                    continue
                pln = ln.strip().split("\t")
                id1, id2 = pln[3], pln[4]
                nodes[id1] = 1
                nodes[id2] = 1
                if id1 not in edges:
                    edges[id1] = {}
                count += 1
                edges[id1][id2] = 1
        print(str(count) + " edges Processed")
        return nodes, edges

    def pruneEqGraph(self, edges):
        uf = UnionFind()
        hash1 = {}
        shash = {}
        num = len(edges)
        count = 0
        eqscores = []
        for u in edges:
            print(count, num, end="\r", flush=True)
            data = []
            scores = []
            for v in edges[u]:
                r = self.rho(u, v, edges, hash1, shash)
                data += [r]
                scores += [[r, v]]
                uf[u], uf[v]
            filter1 = sorted(scores, reverse=True)
            for s in filter1:
                if uf[u] != uf[s[1]]:
                    uf.union(u, s[1])
                    eqscores.append([u, s[0], s[1]])
                    # print(u, s[0], s[1])
                    break
            count += 1
        return pd.DataFrame(eqscores)

    def rho(self, u, v, edges, hash1, shash):
        if shash and u in shash and v in shash[u]:
            return shash[u][v]
        gu = self.gamma(u, edges, hash1)
        gv = self.gamma(v, edges, hash1)
        g_union = set(gu).union(gv)
        g_int = set(gu).intersection(gv)
        res = 0
        if len(g_union) == 0:
            res = 0
        else:
            res = len(g_int) / len(g_union)
        if shash:
            shash[u][v] = res
        return res

    def gamma(self, u, edges, hash1):
        if hash1 and u in hash1:
            return hash1[u]
        res = [u] + list(edges[u].keys())
        if hash1:
            hash1[u] = res
        return res

    def getClusters(self, df, thr=0.5):
        uf = UnionFind()
        edges = {}
        for i in df.index:
            id1 = df[0][i]
            id2 = df[2][i]
            if id1 not in edges:
                edges[id1] = {}
            edges[id1][id1] = df[1][i]
            if id2 not in edges:
                edges[id2] = {}
            edges[id2][id2] = df[1][i]
            uf[id1], uf[id2]
            if df[1][i] > thr:
                uf.union(id1, id2)
                edges[id1][id2] = df[1][i]
                edges[id2][id1] = df[1][i]
        rank = {}
        for k in edges:
            rank[k] = len(edges[k])
        cls = {}
        for k in uf.to_sets():
            l = sorted(k, key=lambda x: rank[x], reverse=True)
            cls[l[0]] = [len(l), l]
        return cls

    def getClustersGraph(self, net, cls):
        nodes = cls
        ids = [k for k in cls]
        eqgraph = []
        count = 0
        for u in ids:
            if (count % 100) == 0:
                print(count, end="\r", flush=True)
            nl = nodes[u][0]
            mid = int(nl / 2)
            l = [u]
            if nl > 2:
                l = [u, nodes[u][1][1], nodes[u][1][mid]]
            if nl > 10:
                l += [nodes[u][1][i] for i in [int(mid / 4), int(mid / 2), mid - 1]]
            ru = self.getNCodes(net, l)
            for c in range(1, 7):
                if c not in ru:
                    continue
                for v in ru[c]:
                    if v not in nodes:
                        continue
                    if nl > 10 and ru[c][v] < 3:
                        continue
                    eqgraph.append([u, str(c), v, ru[c][v]])
            count += 1
        return pd.DataFrame(eqgraph)

    def saveClusters(self, ofile, cls):
        fp = open(ofile, "w")
        for k in sorted(cls, key=lambda x: cls[x][0], reverse=True):
            l1 = [k, cls[k][0]] + cls[k][1]
            l1 = [str(k) for k in l1]
            fp.write("\t".join(l1) + "\n")
        fp.close()
        return

    def getNCodes(net, l):
        relations = {}
        for u in l:
            res = net.readBlockID(u)
            i = net.balancedhash[u]
            for j in range(len(net.balanced)):
                v = net.balanced[j]
                code = res[j]
                if code <= 0:
                    continue
                if code not in relations:
                    relations[code] = {}
                if v not in relations[code]:
                    relations[code][v] = 0
                relations[code][v] += 1
        return relations


def main(file_rl):
    file_base = file_rl.split(".")[0]
    file_res_txt = file_base + "-res.txt"
    _, edges = BIGraph.readEqGraph(file_res_txt)
    df = BIGraph.pruneEqGraph(edges)
    df.to_csv(file_base + "-eq.txt", sep="\t", header=False, index=False)
    cls = BIGraph.getClusters(df)
    BIGraph.saveClusters(file_base + "-cls.txt", cls)
    net = BINetwork(file_rl)
    net.init()
    cg = BIGraph.getClustersGraph(net, cls)
    cg.to_csv(file_base + "-eq-g.txt", sep="\t", header=False, index=False)


if __name__ == "__main__":
    file_rl = sys.argv[1]
    main(file_rl)
