from graph import *
from itertools import permutations, combinations
import numpy as np
import time


# 计算满足特定限制的mapping
def cal_mapping(qubit_list: list, graph: pGraph):
    min_cost = max_dist
    best_mapping = []

    for central_qubit in range(graph.leng):
        neighbors = [i for i in range(graph.leng) if graph.C[central_qubit][i] == 1]
        if len(neighbors) < 3:
            continue
        
        all_conbinations = [list(comb) + [central_qubit] for comb in combinations(neighbors, 3)]
        all_permutations = [list(perm) for i in range(len(all_conbinations)) for perm in permutations(all_conbinations[i])]
        for permutation in all_permutations:
            cost = 0
            for i in range(4):
                cost += graph.C[qubit_list[i]][permutation[i]]
            if cost < min_cost:
                min_cost = cost
                best_mapping = [(qubit_list[i], permutation[i]) for i in range(4)]

    return best_mapping, min_cost


def obstacles_in_path(path: list, in_place: set) -> int:
    cnt = 0
    for i in range(1, len(path)):
        if path[i] in in_place:
            cnt += 1
    return cnt


# 当path中某一步移动到目标结构中后，之后的移动都将在目标结构中进行
def optimize_path(path: list, target_place: set, graph: pGraph):
    found = False
    for i in range(len(path)):
        if not found:
            if path[i] in target_place:
                found = True
        else:
            if path[i] not in target_place:
                for tp in target_place:
                    if graph.C[path[i]][tp] == 1 and graph.C[tp][path[i + 1]] == 1:
                        path[i] = tp
                        return


def is_executable_path(path: list, qubit_bind: dict, qubit_list: list, in_place:set) -> bool:
    for i in range(1, len(path)):
        if qubit_bind[path[i]] in qubit_list and path[i] not in in_place:
            return False
    return True

# 计算满足qubit_list的最优SWAP路径
# 返回目标物理qubit和SWAP路径
# 目标物理qubit为长度为2或4的list

def diag(qubit_list : list, Graph : pGraph, pi : dict):
    C = Graph.C
    adj = Graph.adj
    swaps = []
    cnots = []
    # ses = []
    # if len(pobs) == 2:
    #     path = Graph.get_path(pobs[0], pobs[1])
    #     for i in range(len(path) - 2):
    #         swaps.append((path[i], path[i + 1]))
    #     target_qubits = [path[-1], path[-2]]
    #     return target_qubits, swaps
    def cal_distance(_pi):
        res = 0
        for i in qubit_list:
            for j in qubit_list:
                res += C[_pi[i]][_pi[j]]
        return res
    def get_cnots(targets:list):
        if len(targets) == 2:
            return [[targets[0], targets[1]]]
        v0 = min(targets, key=lambda x: sum([C[x][i] for i in targets]))
        v3 = max(targets, key=lambda x: C[x][v0])
        targets.remove(v0)
        targets.remove(v3)
        v2 = min(targets, key=lambda x: C[x][v3])
        targets.remove(v2)
        v1 = targets[0]
        if v3 in adj[v0]:
            cnots.append([v0, v3])
        else:
            cnots.append([v2, v3])
        cnots.append([v0, v2])
        cnots.append([v0, v1])
    def is_connected():
        pe = [pi[q] for q in qubit_list]
        visited = [pe[0]]
        flag = True
        while(flag):
            flag = False
            for v in visited.copy():
                for _v in adj[v]:
                    if _v in pe and _v not in visited:
                        visited.append(_v)
                        flag = True
        return len(visited) == len(pe)
    while(not is_connected()):
        cand_SWAP = []
        pe = [pi[q] for q in qubit_list]
        for v in pe:
            for _v in adj[v]:
                cand_SWAP.append([v, _v])
        mscore = 1e10
        mpi = {}
        mswap = []
        for swap in cand_SWAP:
            temp_pi = pi.copy()
            for q, v in pi.items():
                if v in swap:
                    temp_pi[q] = swap[1] + swap[0] - v
            score = cal_distance(temp_pi)
            if score < mscore:
                mscore = score
                mpi = temp_pi
                mswap = swap
        pi = mpi
        swaps.append(mswap.copy())
    pe = [pi[q] for q in qubit_list]
    get_cnots(pe.copy())
    return pe, swaps, cnots

def _cal_swap(qubit_list: list, graph: pGraph):
    if not (len(qubit_list) == 2 or len(qubit_list) == 4):
        raise RuntimeError("qubit_list must have 2 or 4 elements")
    
    swaps = []

    if len(qubit_list) == 2:
        path = graph.get_path(qubit_list[0], qubit_list[1])
        for i in range(len(path) - 2):
            swaps.append((path[i], path[i + 1]))
        target_qubits = [path[-1], path[-2]]
        return target_qubits, swaps

    C = graph.C
    vc = -1
    md = 10 ** 10
    for i in range(len(C)):
        d = sum([C[i][q] for q in qubit_list])
        if d < md:
            md = d
            vc = i
    

    return target_qubits, swaps


if __name__ == "__main__":
    start_time = time.time()

    # graph = load_graph("grid0505")
    # qubit_list = [0, 12, 22, 24]

    graph = load_graph("grid0404")
    qubit_list = [0, 1, 8, 10]
    # qubit_list = [0, 125]

    # G = np.array(
    #     [
    #         [0, 1, 0, 0, 0],
    #         [1, 0, 1, 0, 0],
    #         [0, 1, 0, 1, 1],
    #         [0, 0, 1, 0, 0],
    #         [0, 0, 1, 0, 0],
    #     ]
    # )
    # C = np.array(
    #     [
    #         [0, 1, 100000, 100000, 100000],
    #         [1, 0, 1, 100000, 100000],
    #         [100000, 1, 0, 1, 1],
    #         [100000, 100000, 1, 0, 100000],
    #         [100000, 100000, 1, 100000, 0],
    #     ]
    # )
    # dist = floyd_warshall_vectorized(C)
    # graph = pGraph(G, dist)
    # qubit_list = [0, 1, 2, 3]

    graph_time = time.time()

    pi = {i:i for i in range(16)}
    target_qubits, swaps, cnots = diag(qubit_list, graph, pi)
    end_time = time.time()

    print(f"Target qubits: {target_qubits}")
    print(f"Swaps: {swaps}")
    print(f"cnots: {cnots}")
    print(f"Graph time:\t {graph_time - start_time:.2f} s")
    print(f"Calculation time:\t {end_time - graph_time:.2f} s")
    print(f"Total time:\t {end_time - start_time:.2f} s")
