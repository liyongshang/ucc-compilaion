from sys import maxsize
from base import *
from copy import deepcopy
from swap_path import diag
import ast
from tools import showMappingOnGrid, para
import steiner_forest

class TreeNode:
    def __init__(self, number: int, children: list, t, parent=None):
        self.number = number
        self.children = children
        self.t = t
        self.parent = parent
        self.covered = False

    def link(self, child):
        self.children.append(child)
        child.parent = self

    def show(self, layer=0):
        for i in range(layer):
            print("   ", end="")
        print(self.number)
        for c in self.children:
            c.show(layer + 1)

    def traverse(self):
        res = []
        q = [self]
        while q:
            t = q.pop()
            res.append(t)
            q.extend(t.children)
        return res

    def reroot(self):
        if self.parent:
            self.parent.children.remove(self)
            temp = self.parent
            self.parent = None
            temp.reroot()
            self.link(temp)
        return


def findCenter1(graph, V):
    center, minNum = -1, maxsize
    cnt = 0
    for i in V:
        for j in V:
            if i != j:
                cnt += graph[i][j]
        if cnt < minNum:
            center = i
            minNum = cnt
        cnt = 0
    return center


def distToSeleted(dist, V, seleted):
    distance = 0
    for i in V:
        if not seleted[i]:
            minDistance = maxsize
            for j in range(len(dist)):
                if seleted[j] and dist[i][j] < minDistance:
                    minDistance = dist[i][j]
            distance += minDistance
    return distance


def generate(graph, V):
    # 将所有节点转换为TreeNode
    nodes = [TreeNode(n, [], 0) for n in range(len(graph))]
    for i in V:
        nodes[i].covered = True
    vNum = len(graph.data)

    debug = False
    # 寻找中心点
    center = findCenter1(graph.C, V)

    seleted = [False] * vNum
    seleted[center] = True
    cnt = len(V) - 1
    while cnt > 0:
        visited = seleted.copy()
        minCost, parent, child = maxsize, -1, -1
        for i in range(vNum):
            if seleted[i]:
                for temp in graph.data[i].adj:
                    if not visited[temp]:
                        # visited[temp] = True
                        tempSeleted = seleted.copy()
                        tempSeleted[temp] = True
                        distance = distToSeleted(graph.C, V, tempSeleted)
                        cost = (
                            nodes[i].t + distance * (10**3) - int(temp in V) * (10**6)
                        )
                        if cost < minCost:
                            minCost, parent, child = cost, i, temp
        if debug:
            print(parent, "->", child)
            input()
        nodes[parent].link(nodes[child])
        seleted[child] = True
        if child in V:
            cnt -= 1
    return nodes[center]


def generate_tree(graph, V, qubit_list):
    tree = generate(graph, V)
    nodes = tree.traverse()
    min_dist = maxsize

    # 如果qubit_list中有节点在V中，直接返回
    for i in qubit_list:
        for node in nodes:
            if node.number == i:
                node.reroot()
                return node, i

    qubit = -1
    link_node = None
    for i in qubit_list:
        for node in nodes:
            if graph.C[i][node.number] < min_dist:
                qubit = i
                link_node = node
                min_dist = graph.C[i][node.number]
    if qubit == -1:
        raise ValueError("No qubit or link_node")
    link_node.reroot()

    return link_node, qubit


def countNodeNumber(tree: TreeNode):
    res = 0
    for c in tree.children:
        res += countNodeNumber(c)
    res += 1
    return res

class FO:
    def __init__(self, obs):
        obs.sort()
        self.obs = obs
        self.sig = list(range(obs[0] + 1, obs[1]))
        if len(obs) > 2:
            self.sig += list(range(obs[2] + 1, obs[3]))
from mylib import secondCost
from time import time
def unitaryCost(u, pi):
    ans = 0
    pobs = [pi[i] for i in u.obs]
    psig = [pi[i] for i in u.sig]
    t0 = time()
    edges = Graph.SteinerBalanceForest(pobs, psig)
    # print('c1: {}'.format(time() - t0))
    C1 = len(edges) - len(u.sig)
    t0 = time()
    C2 = steiner_forest.secondCost(pobs)
    # real_swaps = diag(pobs, Graph)[1]
    # if C2 != len(real_swaps):
    #     raise ValueError('C2 error: {}!={}'.format(C2, len(real_swaps)))
    # print('c2: {}'.format(time() - t0))
    # input()
    return C1 + C2

def evalPaulis(femi_ops, m, debug=False):
    # showMappingOnGrid(m, 5, 5)
    # print()
    def count_swap(gate_list):
        count = 0
        for gate in gate_list:
            if gate[0] == "SWAP":
                count += 1
        return count
    cost = 0
    for op in femi_ops:
        u = FO(op)
        cost += unitaryCost(u, m)
    return -cost


from math import pi




twoQubitEx = [
    ["RZ", [0], pi / 2],
    ["RY", [1], -pi / 2],
    ["RZ", [1], -pi / 2],
    ["CX", [0, 1]],
    ["RY", [0], para("theta", 0.5)],
    ["RZ", [1], -pi / 2],
    ["CX", [0, 1]],
    ["RY", [0], para("theta", -0.5)],
    ["H", [1]],
    ["CX", [0, 1]],
]

cccRy = [
    ["RY", [0], para("theta", 0.125)],
    ["H", [1]],
    ["CX", [0, 1]],
    ["RY", [0], para("theta", -0.125)],
    ["H", [3]],
    ["CX", [0, 3]],
    ["RY", [0], para("theta", 0.125)],
    ["CX", [0, 1]],
    ["RY", [0], para("theta", -0.125)],
    ["H", [2]],
    ["CX", [0, 2]],
    ["RY", [0], para("theta", 0.125)],
    ["CX", [0, 1]],
    ["RY", [0], para("theta", -0.125)],
    ["CX", [0, 3]],
    ["RY", [0], para("theta", 0.125)],
    ["H", [3]],
    ["CX", [0, 1]],
    ["RY", [0], para("theta", -0.125)],
    ["H", [1]],
    ["CX", [0, 2]],
    ["H", [2]],
]

cccRy_linear = [
    ["RY", [0], para("theta", 0.125)],
    ["H", [1]],
    ["CX", [0, 1]],
    ["RY", [0], para("theta", -0.125)],
    ["H", [3]],
    ["CX", [2, 3]],
    ["CX", [0, 2]],
    ["CX", [2, 3]],
    ["CX", [0, 2]],
    ["RY", [0], para("theta", 0.125)],
    ["CX", [0, 1]],
    ["RY", [0], para("theta", -0.125)],
    ["H", [2]],
    ["CX", [0, 2]],
    ["RY", [0], para("theta", 0.125)],
    ["CX", [0, 1]],
    ["RY", [0], para("theta", -0.125)],
    ["CX", [2, 3]],
    ["CX", [0, 2]],
    ["CX", [2, 3]],
    ["CX", [0, 2]],
    ["RY", [0], para("theta", 0.125)],
    ["H", [3]],
    ["CX", [0, 1]],
    ["RY", [0], para("theta", -0.125)],
    ["H", [1]],
    ["CX", [0, 2]],
    ["H", [2]],
]


def transGate(g, m, sign):
    realG = deepcopy(g)
    realG[1] = [m[i] for i in realG[1]]
    if sign:
        realG[3].coeff *= -1
    return realG


def addCnotTree(graph, electonics: list, V: list, gateList: list):
    return


def addSWAP(graph, electonics: list, gateList: list):
    t = max(electonics, key=lambda e: sum([graph.C[e][i] for i in electonics]))
    electonics.remove(t)
    candidatePos = graph.data[t].adj.copy()


def add_target_qubit(
    graph: pGraph, target_qubit: int, link_node: TreeNode, electonics: list
):
    gates = []
    qubit_bind = {i: i for i in electonics}

    if target_qubit == link_node.number:
        for child in link_node.children:
            if child.number in electonics:
                gates.append(["CZ", [link_node.number, child.number]])
        return gates

    path = graph.get_path(target_qubit, link_node.number)

    for i in range(len(path) - 2):
        gates.append(["SWAP", [path[i], path[i + 1]]])
        if qubit_bind.get(path[i + 1]) != None:
            qubit_bind[path[i]], qubit_bind[path[i + 1]] = (
                qubit_bind[path[i + 1]],
                qubit_bind[path[i]],
            )
        else:
            qubit_bind[path[i + 1]] = qubit_bind[path[i]]
            del qubit_bind[path[i]]
    gates.append(["CZ", [path[-2], link_node.number]])

    electonics.clear()
    for k in qubit_bind:
        electonics.append(k)
    return gates

def eliminateZ(graph: pGraph, electonics: list, V: list, isLeft = True, conAware = True):
    ops = []
    if not conAware:
        ops.append(['CZ', [electonics[-1], V[0]]])
        for i in range(0, len(V) - 1):
            ops.append(['CX', [V[i + 1], V[i]]])
        if isLeft:
            ops.reverse()
        return ops
    edges = Graph.SteinerBalanceForest(electonics, V)
    isSwapped = [False for i in range(nHardwareQubit)]
    for e in edges:
        if e[0] not in electonics + V:
            if isSwapped[e[0]] == False:
                ops.append(['SWAP', list(e)])
                isSwapped[e[0]] = True
            else:
                ops.append(['CX', [e[1], e[0]]])
        if e[0] in electonics:
            ops.append(['CZ', list(e)])
        if e[0] in V:
            ops.append(['CX', [e[1], e[0]]])
    if isLeft:
        ops.reverse()
    return ops
    nodes = [TreeNode(n, [], 0) for n in range(len(graph))]
    for i in V:
        nodes[i].covered = True
    vNum = len(graph.data)

    debug = False
    ops = []

    seleted = [False] * vNum
    for e in electonics:
        seleted[e] = True
    cnt = len(V)
    while cnt > 0:
        visited = seleted.copy()
        minCost, parent, child = maxsize, -1, -1
        for i in range(vNum):
            if seleted[i]:
                for temp in graph.data[i].adj:
                    if not visited[temp]:
                        # visited[temp] = True
                        tempSeleted = seleted.copy()
                        tempSeleted[temp] = True
                        distance = distToSeleted(graph.C, V + electonics, tempSeleted)
                        cost = (
                            nodes[i].t + distance * (10**3) - int(temp in V) * (10**6)
                        )
                        if cost < minCost:
                            minCost, parent, child = cost, i, temp
        if debug:
            print(parent, "->", child)
            # input()
        nodes[parent].link(nodes[child])
        seleted[child] = True
        if child in V:
            cnt -= 1
        if parent in electonics:
            ops.append(['CZ', [child, parent]])
        elif not nodes[parent].covered:
            ops.append(['SWAP', [child, parent]])
            nodes[parent].covered
        else:
            ops.append(['CX', [child, parent]])
    if isLeft:
        ops.reverse()
    return ops




def addDoubleFermionic(graph, op, gateList: list, conAware = True, pi : dict = {}):
    rpi = {i : -1 for i in range(nHardwareQubit)}
    for k, v in pi.items():
        rpi[v] = k

    electonics = op
    V = [i for i in range(op[0] + 1, op[1])]
    if len(op) == 4:
        V += [i for i in range(op[2] + 1, op[3])]
    target_qubits, swap_gates = electonics, []
    if conAware:
        target_qubits, swap_gates, diag_cnots = diag(electonics, graph, pi.copy())
    # _elec = [pi[e] for e in electonics]


    for swap in swap_gates:
        gateList.append(["SWAP", list(swap)])
        Q1, Q2 = swap[0], swap[1]
        q1, q2 = rpi[Q1], rpi[Q2]
        rpi[Q1], rpi[Q2] = q2, q1
        # sg = list(swap)
        # for i in range(len(_elec)):
        #     if _elec[i] == sg[0]:
        #         _elec[i] = sg[1]
        #     elif _elec[i] == sg[1]:
        #         _elec[i] = sg[0]
    temp_pi = {}
    for k,v in rpi.items():
        if v != -1:
            temp_pi[v] = k

    if conAware:
        electonics = [temp_pi[e] for e in electonics]
        V = [temp_pi[k] for k in V]
    if len(V) == 0:
        Zeliminate_gates = []
    else:
        Zeliminate_gates = eliminateZ(graph, electonics, V, True, conAware)

    # if conAware:
    #     electonics = [pi[e] for e in electonics]
    #     V = [pi[k] for k in V]
    # if len(V) == 0:
    #     Zeliminate_gates = []
    # else:
    #     Zeliminate_gates = eliminateZ(graph, electonics, V, True, conAware)
    # 正向添加
    gateList.extend(Zeliminate_gates)

    if len(electonics) == 2:
        signs = [1, 0]
    else:
        signs = [1, 1, 0, 0]
    ms = {}
    for e, s in zip(target_qubits, signs):
        ms[e] = s
    # for i in range(1, len(target_qubits)):
    #     gateList.append(["CX", [target_qubits[0], target_qubits[i]]])
    #     c = target_qubits[0]
    #     t = target_qubits[i]
    #     if not ms[c] ^ ms[t]:
    #         gateList.append(['X', [t]])
    for cnot in diag_cnots:
        c, t = cnot[0], cnot[1]
        gateList.append(["CX", [c, t]])
        if not ms[c] ^ ms[t]:
            gateList.append(['X', [t]])

    # 添加受控RY
    if len(target_qubits) == 2:
        for g in twoQubitEx:
            gateList.append(transGate(g, target_qubits, False))
    elif len(target_qubits) == 4:
        target_qubits = [cnot[1] for cnot in diag_cnots] + [diag_cnots[-1][0]]
        target_qubits.reverse()
        if diag_cnots[0][0] == diag_cnots[-1][0]:
            for g in cccRy:
                gateList.append(transGate(g, target_qubits, False))
        else:
            for g in cccRy_linear:
                gateList.append(transGate(g, target_qubits, False))
    else:
        raise ValueError("The number of target qubits is not 2 or 4")

    # 反向添加
    # for i in range(len(target_qubits) - 1, 0, -1):
    #     c = target_qubits[0]
    #     t = target_qubits[i]
    #     if not ms[c] ^ ms[t]:
    #         gateList.append(['X', [t]])
    #     gateList.append(["CX", [target_qubits[0], target_qubits[i]]])
    diag_cnots.reverse()
    for cnot in diag_cnots:
        c, t = cnot[0], cnot[1]
        if not ms[c] ^ ms[t]:
            gateList.append(['X', [t]])
        gateList.append(["CX", [c, t]])
        
    for i in range(len(Zeliminate_gates) - 1, -1, -1):
        gateList.append(Zeliminate_gates[i])
    for i in range(len(swap_gates) - 1, -1, -1):
        gateList.append(["SWAP", list(swap_gates[i])])

    return


if __name__ == "__main__":
    print(1)
    start_time = time()

    graph = load_graph("grid0404")
    steiner_forest.set_graph_data(Graph.adj, Graph.C, nHardwareQubit)
    qubit_list = [0, 1, 8, 10]

    graph_time = time()

    pi = {i:i for i in range(16)}
    gateList = []
    addDoubleFermionic(graph, qubit_list, [9], gateList, True, pi)
    end_time = time()

    for g in gateList:
        print(g)
    print(f"Graph time:\t {graph_time - start_time:.2f} s")
    print(f"Calculation time:\t {end_time - graph_time:.2f} s")
    print(f"Total time:\t {end_time - start_time:.2f} s")