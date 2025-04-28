__all__ = [
    "load_graph",
    "load_coupling_map",
    "pNode",
    "pGraph",
    "dijkstra",
    "max_dist",
    "graph_from_coupling",
    "floyd_warshall_vectorized",
]
import os
import numpy as np
import sys
import csv
import json
from copy import deepcopy
import random
import steiner_forest

max_dist = 100000
max_size = 100000
package_directory = os.path.dirname(os.path.abspath(__file__))


class pNode:
    def __init__(self, idx):
        # self.child = []
        self.idx = idx
        self.adj = []
        self.lqb = None  # logical qubit
        # self.parent = []

    def add_adjacent(self, idx):
        self.adj.append(idx)


class pGraph:
    def __init__(self, G, C):
        n = G.shape[0]
        self.leng = n
        self.G = G  # adj matrix
        self.C = C  # cost matrix
        self.data = []
        self.coupling_map = []
        self.adj = [[] for i in range(n)]
        for i in range(n):
            nd = pNode(i)
            for j in range(n):
                if G[i, j] == 1:
                    nd.add_adjacent(j)
                    self.coupling_map.append([i, j])
                    self.adj[i].append(j)
                    # self.adj[j].append(i)
            self.data.append(nd)

    def find_center(self):
        c = -1
        md = 10 ** 10
        for i in range(len(self.C)):
            d = max(self.C[i])
            if d < md:
                md = d
                c = i
        return c

    def __getitem__(self, idx):
        return self.data[idx]

    def __len__(self):
        return self.leng

    def copy(self):
        pgh = pGraph(self.G, self.C)
        for i in range(len(self.data)):
            pgh.data[i].lqb = self.data[i].lqb
        return pgh

    def print(self):
        for row in self.G:
            print(row)

    def get_path(self, start, end):
        dist = self.C
        if dist[start][end] == max_dist:
            raise RuntimeError("No path between qubits")
        path = [start]
        while start != end:
            for i in range(self.leng):
                if dist[start][i] == 1 and dist[i][end] == dist[start][end] - 1:
                    path.append(i)
                    start = i
                    break
        return path

    def SteinerBalanceForest(self, roots : list, A : list):
        edges = steiner_forest.SteinerBalanceForest(roots, A)
        return edges
        

def minDistance(dist, sptSet):
    minv = max_size
    min_index = -1
    n = len(dist)
    for v in range(n):
        if dist[v] < minv and sptSet[v] == False:
            minv = dist[v]
            min_index = v
    return min_index


def dijkstra(dist_matrix, src):
    n = dist_matrix.shape[0]
    dist = [dist_matrix[src, i] for i in range(n)]
    sptSet = [False] * n
    for cout in range(n):
        u = minDistance(dist, sptSet)
        if u == -1:
            break
        sptSet[u] = True
        for v in range(n):
            if (
                dist_matrix[u, v] > 0
                and sptSet[v] == False
                and dist[v] > dist[u] + dist_matrix[u, v]
            ):
                dist[v] = dist[u] + dist_matrix[u, v]
    for i in range(n):
        dist_matrix[src, i] = dist[i]
        dist_matrix[i, src] = dist[i]


def check_and_convert_adjacency_matrix(dist_matrix):
    matrix = np.asarray(dist_matrix)

    nrows, ncols = matrix.shape
    assert nrows == ncols
    n = nrows

    assert (np.diagonal(matrix) == 0).all()

    return matrix, n


def floyd_warshall_vectorized(dist_matrix):
    matrix, n = check_and_convert_adjacency_matrix(dist_matrix)

    for k in range(n):
        matrix = np.minimum(matrix, matrix[np.newaxis, k, :] + matrix[:, k, np.newaxis])

    return matrix


def is_code_reduced(code):
    if code in ["melbourne", "mahattan"]:
        reduced = True
    else:
        reduced = False
    return reduced


def load_graph(code, dist_comp=True, len_func=lambda x: x):
    if "grid" in code:
        a = int(code[4:6])  # a*b的阵列
        b = int(code[6:8])
        n = a * b
        G = np.zeros((n, n), dtype=int)
        C = np.ones((n, n), dtype=int) * max_dist
        for i in range(a):
            for j in range(b):
                x = i * b + j
                C[x][x] = 0
                if j < b - 1:
                    G[x][x + 1] = G[x + 1][x] = 1
                    C[x][x + 1] = C[x + 1][x] = 1
                if i < a - 1:
                    G[x][x + b] = G[x + b][x] = 1
                    C[x][x + b] = C[x + b][x] = 1
        # G[10][11] = G[11][10] = 0
        # C[10][11] = C[11][10] = max_dist
        if dist_comp == True:
            dist = floyd_warshall_vectorized(C)
        return pGraph(G, dist)
    if code == "sycamore":
        adj = json.load(open(os.path.join(package_directory, "sycamore.json")))[
            "adj"
        ]
        n = len(adj)
        G = np.zeros((n, n), dtype=int)
        C = np.ones((n, n), dtype=int) * max_dist
        for i in range(n):
            C[i][i] = 0
        for i in range(n):
            for a in adj[i]:
                j = int(a["v"][2:-1])
                G[i][j] = G[j][i] = 1
                C[i][j] = C[j][i] = 1
            #     print(j, ' ', end="")
            # print('\n')
        if dist_comp == True:
            dist = floyd_warshall_vectorized(C)
        return pGraph(G, dist)

    reduced = is_code_reduced(code)
    pth = os.path.join(package_directory, "data", "ibmq_" + code + "_calibrations.csv")
    cgs = []
    n = 0
    with open(pth, "r") as cf:
        g = csv.DictReader(cf, delimiter=",", quotechar='"')
        for i in g:
            cxval = ""
            for j in i.keys():
                if j.find("CNOT") != -1:
                    cxval = i[j]
            n += 1
            if ";" in cxval:
                dc = ";"
            else:
                dc = ","
            for j in cxval.split(dc):
                cgs.append(j.strip())  # 所有CX及其错误率
    if reduced:
        n -= 1  # qubit n is not used
    G = np.zeros((n, n))
    C = np.ones((n, n)) * max_dist
    for i in range(n):
        C[i, i] = 0
    for i in cgs:
        si1 = i.find("_")
        si2 = i.find(":")
        offset = 0
        if i[:2] == "cx":
            offset = 2
        iq1 = int(i[offset:si1])
        iq2 = int(i[si1 + 1 : si2])
        acc = float(i[si2 + 2 :]) * 1000
        if (iq1 < n and iq2 < n) or not reduced:
            G[iq1, iq2] = 1
            C[iq1, iq2] = 1  # len_func(acc/1000) #q1, q2错误率
    if dist_comp == True:
        dist = floyd_warshall_vectorized(C)
    return pGraph(G, dist)


def graph_from_coupling(coup, dist_comp=True):
    n = max([max(i) for i in coup]) + 1
    G = np.zeros((n, n))
    C = np.ones((n, n)) * max_dist
    for i in range(n):
        C[i, i] = 0
    for i in coup:
        G[i[0], i[1]] = 1
        C[i[0], i[1]] = 1
    if dist_comp == True:
        dist = floyd_warshall_vectorized(C)
    return pGraph(G, dist)


def load_coupling_map(code):
    reduced = is_code_reduced(code)
    pth = os.path.join(package_directory, "data", "ibmq_" + code + "_calibrations.csv")
    cgs = []
    n = 0
    with open(pth, "r") as cf:
        g = csv.DictReader(cf, delimiter=",", quotechar='"')
        for i in g:
            cxval = ""
            for j in i.keys():
                if j.find("CNOT") != -1:
                    cxval = i[j]
            n += 1
            if ";" in cxval:
                dc = ";"
            else:
                dc = ","
            for j in cxval.split(dc):
                cgs.append(j.strip())
    coupling = []
    if reduced:
        n -= 1
    for i in cgs:
        si1 = i.find("_")
        si2 = i.find(":")
        offset = 0
        if i[:2] == "cx":
            offset = 2
        iq1 = int(i[offset:si1])
        iq2 = int(i[si1 + 1 : si2])
        if (iq1 < n and iq2 < n) or not reduced:
            coupling.append([iq1, iq2])
    return coupling
