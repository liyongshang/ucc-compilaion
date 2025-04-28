from qiskit import QuantumCircuit
import numpy as np
import sys

max_size = 10**9
def dummy_mapping(n):
    pi = {}
    for i in range(n):
        pi[i] = i
    return pi


def print_mapping(mapping, ac):
    if 'grid' in ac:
        a = int(ac[4:6])  # a*b的阵列
        b = int(ac[6:8])
        for i in range(a):
            for j in range(b):
                if i % 2 == 1:
                    j = b - 1 - j
                print(mapping[i * b + j], ' ', end="")
            print('\n')
    else:
        rows = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, -1],
                [10, -1, -1, -1, 11, -1, -1, -1, 12, -1, -1],
                list(range(13, 24)),
                [-1, -1, 24, -1, -1, -1, 25, -1, -1, -1, 26],
                list(range(27, 38)),
                [38, -1, -1, -1, 39, -1, -1, -1, 40, -1, -1],
                list(range(41, 52)),
                [-1, -1, 52, -1, -1, -1, 53, -1, -1, -1, 54],
                [-1] + list(range(55, 65))]
        for r in rows:
            for a in r:
                if a != -1:
                    if mapping[a] != -1:
                        print(mapping[a], ' ', end="")
                    else:
                        print('* ', end="")
                else:
                    print(' ', ' ', end="")
            print('\n')

def ps2nodes(ps):
    r = []
    for i in range(len(ps)):
        if ps[i] != 'I':
            r.append(i)
    return r

def compute_block_cover(pauli_block):
    blo_cov = []
    for i in pauli_block:
        for l in ps2nodes(i.ps):
            if l not in blo_cov:
                blo_cov.append(l)
    return blo_cov

def find_path(graph, pid0, pid1, ins):
    minid = -1
    mindist = max_size
    for i in graph[pid0].adj:
        if graph.C[i, pid1] < mindist:
            minid = i
            mindist = graph.C[i, pid1]
    #print(pid0, pid1,minid)
    if minid == pid1:
        ins.append(['c', pid0, minid])
    else:
        ins.append(['s', pid0, minid])
        find_path(graph, minid, pid1, ins)
class TreeNode:
    def __init__(self, number, children:list, t):
        self.number = number
        self.children = children
        self.t = t
        self.covered = False

    def link(self, child):
        self.children.append(child)

    def unlink(self, child):
        self.children.remove(child)

def print_tree1(tree:TreeNode):
    for child in tree.children:
        print(tree.number, '->', child.number)
    for child in tree.children:
        print_tree1(child)

def print_tree2(tree):
    for child in tree.childs:
        print(tree.pid, '->', child.pid)
    for child in tree.childs:
        print_tree2(child)

# 以最小平均距离为标准寻找中心点
def findCenter1(graph, P):
    center, minNum = -1, sys.maxsize
    cnt = 0
    for i in P:
        for j in P:
            if i != j:
                cnt += graph[i][j]
        if cnt < minNum:
            center = i
            minNum = cnt
        cnt = 0
    return center

# 以最小最远距离为标准寻找中心点
def findCenter2(graph, P):
    maxNum, center = sys.maxsize, -1
    for i in P:
        maxTemp = -1
        for j in P:
            if i != j:
                if maxTemp < graph[i][j]:
                    maxTemp = graph[i][j]
        if maxNum > maxTemp:
            maxNum = maxTemp
            center = i
    return center


def distToSeleted(dist, P, seleted):
    distance = 0
    for i in P:
        if not seleted[i]:
            minDistance = max_size
            for j in range(len(dist)):
                if seleted[j] and dist[i][j] < minDistance:
                    minDistance = dist[i][j]
            distance += minDistance
    return distance

def generate(graph, P):
    # 将所有节点转换为TreeNode
    nodes = [TreeNode(n, [], 0) for n in range(len(graph))]
    for i in P:
        nodes[i].covered = True
    vNum = len(graph.data)

    debug = False
    # 寻找中心点
    center = findCenter1(graph.C, P)

    seleted = [False] * vNum
    seleted[center] = True
    cnt = len(P) - 1
    while cnt > 0:
        visited = seleted.copy()
        minCost, parent, child = max_size, -1, -1
        for i in range(vNum):
            if seleted[i]:
                for temp in graph.data[i].adj:
                    if not visited[temp]:
                        # visited[temp] = True
                        tempSeleted = seleted.copy()
                        tempSeleted[temp] = True
                        distance = distToSeleted(graph.C, P, tempSeleted)
                        cost = nodes[i].t + distance * (10**3) - int(temp in P) * (10**6)
                        if (cost < minCost):
                            minCost, parent, child = cost, i, temp
        if debug:
            print(parent, '->', child)
            input()
        nodes[parent].link(nodes[child])
        seleted[child] = True
        if child in P:
            cnt -= 1
    return nodes[center]

def cal_depth(root:TreeNode):
    root.t = 0
    if len(root.children) == 0:
        return
    for child in root.children:
        cal_depth(child)
    root.children.sort(key=lambda x:x.t)
    for i in range(len(root.children)):
        root.t = max(root.t, root.children[i].t) + 1
    return

def generate1(graph, P):
    # 将所有节点转换为TreeNode
    nodes = [TreeNode(n, [], 0) for n in range(len(graph))]
    for i in P:
        nodes[i].covered = True
    vNum = len(graph.data)

    debug = False
    # 寻找中心点
    center = findCenter1(graph.C, P)

    seleted = [False] * vNum
    seleted[center] = True
    cnt = len(P) - 1
    while cnt > 0:
        visited = seleted.copy()
        minCost, parent, child = max_size, -1, -1
        for i in range(vNum):
            if seleted[i]:
                for temp in graph.data[i].adj:
                    if not visited[temp]:
                        # visited[temp] = True
                        tempSeleted = seleted.copy()
                        tempSeleted[temp] = True
                        distance = distToSeleted(graph.C, P, tempSeleted)
                        nodes[i].link(nodes[temp])
                        cal_depth(nodes[center])
                        cost = nodes[center].t + distance * (10**3) - int(temp in P) * (10**6)
                        # if P == [2,1,4,5,3,7]:
                        #     print('({},{})'.format(i, temp), ' : ', nodes[center].t, ', ', cost)
                        nodes[i].unlink(nodes[temp])
                        if (cost < minCost):
                            minCost, parent, child = cost, i, temp
        if debug:
            print(parent, '->', child)
            input()
        nodes[parent].link(nodes[child])
        seleted[child] = True
        if P == [2,1,4,5,3,7]:
            print(seleted)
            print(parent, '->', child)
        if child in P:
            cnt -= 1
    print('my tree:')
    print_tree1(nodes[center])
    return nodes[center]

def go_synthesis(graph, tree, ins):
    if len(tree.children) > 0:
        for i in range(len(tree.children)):
            go_synthesis(graph, tree.children[i], ins)
            # tree.children[i].t += graph.C[tree.number][tree.children[i].number]
        tree.children.sort(key=lambda x: x.t)
        tree.t = max([tree.children[i].t + i + 1 for i in range(len(tree.children))])
        for child in tree.children:
            # if graph.G[tree.number][child.number] == 0:
            #     find_path(graph, child.number, tree.number, ins)
            # else:
            #     ins.append(['c', child.number, tree.number])
            if tree.covered:
                ins.append(['c', child.number, tree.number])
            else:
                ins.append(['s', child.number, tree.number])
                tree.covered = True
    else:
        tree.t = 0

def go_synthesis1(graph, qc, psd, pi, param):
    res = 0
    nodes = [pi[node] for node in ps2nodes(psd)]
    tree = generate(graph, nodes)
    for i in range(len(psd)):
        if psd[i] == 'X':
            # qc.u(np.pi / 2, 0, np.pi, pi[i])
            qc.h(pi[i])
        elif psd[i] == 'Y':
            # qc.u(np.pi / 2, -np.pi / 2, np.pi / 2, pi[i])
            qc.rx(np.pi/2, pi[i])
    ins = []
    go_synthesis(graph, tree, ins)
    for ins1 in ins:
        if ins1[0] == 's':
            res += 1
            qc.swap(ins1[1], ins1[2])
        else:
            qc.cx(ins1[1], ins1[2])
    qc.rz(param, tree.number)
    ins.reverse()
    for ins1 in ins:
        if ins1[0] == 's':
            res += 1
            qc.swap(ins1[1], ins1[2])
        else:
            qc.cx(ins1[1], ins1[2])
    for i in range(len(psd)):
        if psd[i] == 'X':
            # qc.u(np.pi / 2, 0, np.pi, pi[i])
            qc.h(pi[i])
        elif psd[i] == 'Y':
            # qc.u(np.pi / 2, -np.pi / 2, np.pi / 2, pi[i])
            qc.rx(-np.pi/2, pi[i])
    return res