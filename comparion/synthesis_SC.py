from synthesis_sd import *
from arch import *
#from synthesis_FT import assign_time_parameter
from functools import partial
from functions import go_synthesis1, print_tree2

def assign_time_parameter(ps_layers, time_parameter):
    for i in ps_layers:
        for j in i:
            for k in range(len(j)):
                j[k].real += time_parameter
                j[k].coeff += time_parameter

def dummy_local_move(qc, graph, pauli_map, src, target):
    an = graph[src]
    ml = 10000
    mn = -1
    for i in an.adj:
        if graph.C[i, target] < ml:
            ml = graph.C[i, target]
            mn = i
    if ml == 0: # we find it
        qc.cx(src, target)
    else:
        swap_nodes(pauli_map, graph[src], graph[mn])
        qc.swap(src, mn)
        dummy_local_move(qc, graph, pauli_map, mn, target)

from qiskit import QuantumCircuit

class treeNode:
    def __init__(self, pid, status):
        self.pid = pid
        self.status = 0 # 0 means inner node; 1 means leaf node

class tree: #深度优先树
    def __init__(self, graph, dp, parent=None, depth=0):
        self.childs = []
        self.leaf = []
        self.depth = depth
        self.pid = dp[0] #根节点在graph中的编号
        self.parent = self
        if len(dp) == 1:
            self.status = 1
            self.leaf = [self]
        else:
            self.status = 0   
            st = []
            # for i in range(1, len(dp)):
            #     if dp[i] in graph[self.pid].adj:
            #         st.append(i)
            # st.append(len(dp))
            # for i in range(len(st)-1):
            #     child = tree(graph, dp[st[i]:st[i+1]], parent=self, depth=self.depth+1)
            #     self.childs.append(child)
            #     self.leaf += child.leaf
            st.append(1)
            p = 2
            while p < len(dp):
                f = False
                for i in range(st[len(st) - 1], p):
                    if dp[p] in graph[dp[i]].adj:
                        f = True
                        break
                if f == False:
                    st.append(p)
                p += 1
            st.append(len(dp))
            for i in range(len(st)-1):
                child = tree(graph, dp[st[i]:st[i+1]], parent=self, depth=self.depth+1)
                self.childs.append(child)
                self.leaf += child.leaf
        if parent != None:
            self.parent = parent

# swap tree nodes only change its physical qubit id and logical id mapping;
def swap_tree_node(t0, t1):
    pass

def pauli_single_gates(qc, pauli_map, ps, left=True):
    if left == True:
        for i in range(len(ps)):
            if ps[i] == 'X':
                # qc.u(np.pi/2, 0, np.pi, pauli_map[i])
                qc.h(pauli_map[i])
            elif ps[i] == 'Y':
                # qc.u(np.pi/2, -np.pi/2, np.pi/2, pauli_map[i])
                qc.rx(np.pi/2, pauli_map[i])
    else:
        for i in range(len(ps)):
            if ps[i] == 'X':
                # qc.u(np.pi/2, 0, np.pi, pauli_map[i])
                qc.h(pauli_map[i])
            elif ps[i] == 'Y':
                # qc.u(-np.pi/2, -np.pi/2, np.pi/2, pauli_map[i])
                qc.rx(-np.pi/2, pauli_map[i])

def tree_synthesis1(qc, graph, pauli_map, ptree, psd, param): # psd : pauli string
    ps = psd.ps
    psn = ps2nodes(ps)
    pauli_single_gates(qc, pauli_map, ps, left=True) # add single gates
    lfs = ptree.leaf # first in, first out
    swaps = {}
    cnum = len(psn) - 1
    lc = 0 # number of CX
    cxs = []
    # 自底向上，生成左边
    while lfs != []:
        lfs = sorted(lfs, key=lambda x: -x.depth)
        l = lfs[0]
        #print(l.parent.pid)
        if l.depth == 0:
            # psd.real may be zero
            qc.rz(param, l.pid) # qc.rz(2*psd.real+1, l.pid)
            break
        # actually, if psn is empty in the middle, we can stop it first
        # and the choice of root is also important
        if graph[l.pid].lqb in psn:
            if graph[l.parent.pid].lqb in psn:
                qc.cx(l.pid, l.parent.pid)
                c1 = []
                c1.append(l.pid)
                c1.append(l.parent.pid)
                cxs.append(c1)
                lc += 1
            else:
                qc.swap(l.pid, l.parent.pid)
                swaps[l.parent.pid] = l
                swap_nodes(pauli_map, graph[l.pid], graph[l.parent.pid])
        else:
            pass #lfs.remove(l)
        if l.parent not in lfs:
                lfs.append(l.parent)
        lfs.remove(l)
        # print(lfs)
    if lc != cnum:
        print('lala left:',psd.ps, cnum, lc)
        pqs = [pauli_map[i] for i in psn]
        print(pqs)
        print('tree:')
        lfs = ptree.leaf
        while lfs != []:
            lfs = sorted(lfs, key=lambda x: -x.depth)
            l = lfs[0]
            print(l.pid, ' ', l.parent.pid)
            if l.parent not in lfs:
                lfs.append(l.parent)
            lfs.remove(l)
        print('cx:')
        for c in cxs:
            print(c[0], c[1])
        input()
    # 自顶向下，生成右边
    lfs = [ptree]
    rc = 0
    while lfs != []:
        l = lfs[0]
        for i in l.childs:
            if graph[i.pid].lqb in psn:
                qc.cx(i.pid, l.pid)
                rc += 1
                lfs.append(i)
        if l.pid in swaps.keys():
            qc.swap(l.pid, swaps[l.pid].pid)
            swap_nodes(pauli_map, graph[l.pid], graph[swaps[l.pid].pid])
            lfs.append(swaps[l.pid])
        lfs = lfs[1:]
    if rc != cnum:
        print('lala left:',psd.ps, cnum, rc)
    pauli_single_gates(qc, pauli_map, ps, left=False)
    return qc, 2 * len(swaps)

def synthesis_initial(pauli_layers, pauli_map=None, graph=None, qc=None, arch='manhattan'):
    # assign_time_parameter(pauli_layers, 1) #配置\delta t
    lnq = len(pauli_layers[0][0][0]) # logical qubits
    if graph == None:
        G, C = load_graph(arch, dist_comp=True) # G is adj, C is dist
        graph = pGraph(G, C)
    if pauli_map == None:
        pauli_map = dummy_qubit_mapping(graph, lnq)
    else:
        add_pauli_map(graph, pauli_map)
    pnq = len(graph) # physical qubits
    if qc == None:
        qc = QuantumCircuit(pnq)
    return pauli_map, graph, qc

def inter_synthesis(pauli_layers, pauli_map=None, graph=None, qc=None, arch='manhattan'):
    pauli_map, graph, qc = synthesis_initial(pauli_layers, pauli_map, graph, qc, arch)
    for i1 in pauli_layers:
        for i2 in i1:
            pass
from qiskit.circuit import Parameter
def block_opt_SC(pauli_layers, pauli_map=None, graph=None, qc=None, arch='manhattan', synthesis_opt = False):
    pauli_map, graph, qc = synthesis_initial(pauli_layers, pauli_map, graph, qc, arch)
    remain_layers = []
    # print(pauli_layers)
    dp = []
    ins = []
    # lys code
    inner_swaps = 0
    outer_swaps = 0
    params_num = 0
    # lys code
    for i1 in pauli_layers: # i1 : a layer
        for i2 in i1[:1]: # i2 : a block
            # small blocks postponed
            if max([len(ps2nodes(i3.ps)) for i3 in i2]) < 3: # i3 : a string; ps2nodes() -> number of non I
                remain_layers.append([i2])
                continue
            lcover = compute_block_cover(i2) # i2 string 的并集
            itir = compute_block_interior(i2) # i2 string 的交集
            pcover = logical_list_physical(pauli_map, lcover)
            ptir = logical_list_physical(pauli_map, itir)
            # 遍历ptir，找到最大生成树
            lmc = -1 # 生成树的节点数
            lmi = -1 # 生成树的根（物理比特）的编号
            lmt = [] # 生成树的DFS序
            for i3 in ptir: # pauli_map, pcover
                dp = max_dfs_tree(graph, pcover, graph[i3])
                if len(dp) > lmc:
                    lmc = len(dp)
                    lmi = i3
                    lmt = dp
            if len(lmt) == 0:
                lmt = [pcover[0]]
            lcover1 = physical_list_logical(graph, lmt) # 生成树的逻辑比特
            nc = [] # 不在生成树的逻辑比特
            for i3 in lcover:
                if i3 not in lcover1:
                    nc.append(i3)
            # print(nc)
            # 将不连通的逻辑比特移动到和生成树连接
            ins = []
            while nc != []:
                id0, id1 = find_short_node(graph, pauli_map, nc, lmt)
                # print(id0, id1)
                connect_node(graph, pauli_map, pauli_map[id0], pauli_map[id1], ins)
                # do we need to update lmt?
                lmt.append(pauli_map[id0])
                nc.remove(id0)
            for i3 in ins:
                if i3[0] == 'swap':
                    qc.swap(i3[1][0], i3[1][1])
            pcover = logical_list_physical(pauli_map, lcover) # 连接后block涉及的物理比特
            # root is lmi
            dp = max_dfs_tree(graph, pcover, graph[lmi])
            # for i3 in pcover:
            #     if i3 not in dp:
            #         print('e1')
            #         input()
            dt = tree(graph, dp)
            # tns = []
            # lfs = dt.leaf
            # while lfs != []:
            #     lfs = sorted(lfs, key=lambda x: -x.depth)
            #     l = lfs[0]
            #     #print(l.pid, ' ', l.parent.pid)
            #     if l.parent not in lfs:
            #         lfs.append(l.parent)
            #     lfs.remove(l)
            #     tns.append(l.pid)
            # for i3 in dp:
            #     if i3 not in tns:
            #         print(tns)
            #         print(dp)
            #         print('e2')
            #         input()
            parameter = Parameter('phi' + str(params_num))
            for i3 in i2:
                if synthesis_opt:
                    ns = go_synthesis1(graph, qc, i3.ps, pauli_map, i3.coeff * parameter)
                    inner_swaps += ns
                else:
                    _qc, ns = tree_synthesis1(qc, graph, pauli_map, dt, i3, i3.coeff * parameter)
                    inner_swaps += ns
            params_num += 1
        xlist = dp # 第一个block涉及的物理比特
        move_overhead = len(ins) # 外部SWAP个数
        outer_swaps += move_overhead
        for i2 in i1[1:]:
            if max([len(ps2nodes(i3.ps)) for i3 in i2]) < 3:
                remain_layers.append([i2])
                continue
            lcover = compute_block_cover(i2)
            itir = compute_block_interior(i2)
            pcover = logical_list_physical(pauli_map, lcover)
            ptir = logical_list_physical(pauli_map, itir)
            lmc = -1
            lmi = -1
            lmt = []
            for i3 in ptir:
                dp = max_dfs_tree(graph, pcover, graph[i3])
                if len(dp) > lmc:
                    lmc = len(dp)
                    lmi = i3
                    lmt = dp
            if len(lmt) == 0:
                lmt = [pcover[0]]
            lcover1 = physical_list_logical(graph, lmt)
            nc = [] # 除生成树外的逻辑比特
            for i3 in lcover:
                if i3 not in lcover1:
                    nc.append(i3)
            ret = 0
            ins_try = []
            nc_try = nc.copy()
            lmt_try = lmt.copy()
            graph_try = graph.copy()
            pauli_map_try = pauli_map.copy()
            while nc_try != []:
                id0, id1 = find_short_node(graph_try, pauli_map_try, nc_try, lmt_try)
                # 如果路由路径经过第一个block的区域，则ret = -1
                ret = try_connect_node_2(graph_try, pauli_map_try, pauli_map_try[id0], pauli_map_try[id1], ins_try, xlist)
                if ret == -1:
                    remain_layers.append([i2])
                    break
                lmt_try.append(pauli_map_try[id0])
                nc_try.remove(id0)
            if ret == -1:
                continue
            if len(ins_try) > move_overhead:
                remain_layers.append([i2])
                continue
            ins = []
            while nc != []:
                id0, id1 = find_short_node(graph, pauli_map, nc, lmt)
                connect_node(graph, pauli_map, pauli_map[id0], pauli_map[id1], ins)
                lmt.append(pauli_map[id0])
                nc.remove(id0)
            for i3 in ins:
                if i3[0] == 'swap':
                    qc.swap(i3[1][0], i3[1][1])
            outer_swaps += len(ins)
            pcover = logical_list_physical(pauli_map, lcover)
            dp = max_dfs_tree(graph, pcover, graph[lmi])
            dt = tree(graph, dp)
            parameter = Parameter('phi' + str(params_num))
            for i3 in i2:
                if synthesis_opt:
                    ns = go_synthesis1(graph, qc, i3.ps, pauli_map, i3.coeff * parameter)
                    inner_swaps += ns
                else:
                    _qc, ns = tree_synthesis1(qc, graph, pauli_map, dt, i3, i3.coeff * parameter)
                    inner_swaps += ns
            params_num += 1
    # print(remain_layers)
    if remain_layers != []:
        def __key(cost_matrix, pauli_map, ly):
            ns = ps2nodes(ly[0][0].ps)
            ns_len = len(ns)
            if ns_len == 1:
                return 0
            s = 0
            for i in range(ns_len):
                for j in range(i+1,ns_len):
                    s += cost_matrix[pauli_map[ns[i]], pauli_map[ns[j]]]
            return s
        while remain_layers != []:
            # print(remain_layers)
            remain_layers = sorted(remain_layers, key=partial(__key, graph.C, pauli_map))
            picked = remain_layers[0]
            remain_layers = remain_layers[1:]
            for i2 in picked:
                lcover = compute_block_cover(i2)
                itir = compute_block_interior(i2)
                pcover = logical_list_physical(pauli_map, lcover)
                ptir = logical_list_physical(pauli_map, itir)
                lmc = -1
                lmi = -1
                lmt = []
                for i3 in ptir: # pauli_map, pcover
                    dp = max_dfs_tree(graph, pcover, graph[i3])
                    if len(dp) > lmc:
                        lmc = len(dp)
                        lmi = i3
                        lmt = dp
                if len(lmt) == 0:
                    lmt = [pcover[0]]
                lcover1 = physical_list_logical(graph, lmt)
                nc = []
                for i3 in lcover:
                    if i3 not in lcover1:
                        nc.append(i3)
                # print(nc)
                ins = []
                while nc != []:
                    id0, id1 = find_short_node(graph, pauli_map, nc, lmt)
                    # print(id0, id1)
                    connect_node(graph, pauli_map, pauli_map[id0], pauli_map[id1], ins)
                    # do we need to update lmt?
                    lmt.append(pauli_map[id0])
                    nc.remove(id0)
                for i3 in ins:
                    if i3[0] == 'swap':
                        qc.swap(i3[1][0], i3[1][1])
                outer_swaps += len(ins)
                pcover = logical_list_physical(pauli_map, lcover)
                # root is lmi
                dp = max_dfs_tree(graph, pcover, graph[lmi])
                # print(dp)
                dt = tree(graph, dp)
                parameter = Parameter('phi' + str(params_num))
                for i3 in i2:
                    if synthesis_opt:
                        ns = go_synthesis1(graph, qc, i3.ps, pauli_map, i3.coeff * parameter)
                        inner_swaps += ns
                    else:
                        _qc, ns = tree_synthesis1(qc, graph, pauli_map, dt, i3, i3.coeff * parameter)
                        inner_swaps += ns
                params_num += 1
    # print('inner swaps: ', inner_swaps)
    # print('outer_swaps: ', outer_swaps)
    return qc, inner_swaps, outer_swaps

def connected_tree_synthesis(pauli_layers, pauli_map=None, graph=None, qc=None, arch='manhattan'):
    lnq = len(pauli_layers[0][0][0]) # logical qubits
    if graph == None:
        G, C = load_graph(arch, dist_comp=True) # G is adj, C is dist
        graph = pGraph(G, C)
    if pauli_map == None:
        pauli_map = dummy_qubit_mapping(graph, lnq)
    else:
        add_pauli_map(graph, pauli_map)
    pnq = len(graph) # physical qubits
    if qc == None:
        qc = QuantumCircuit(pnq)
    for i1 in pauli_layers:
        for i2 in i1:
            lcover = compute_block_cover(i2)
            pcover = logical_list_physical(pauli_map, lcover)
            lmc = -1
            lmi = -1
            lmt = []
            for i3 in pauli_map:
                dp = max_dfs_tree(graph, pcover, graph[i3])
                if len(dp) > lmc:
                    lmc = len(dp)
                    lmi = i3
                    lmt = dp
            lcover1 = physical_list_logical(graph, lmt)
            nc = []
            for i3 in lcover:
                if i3 not in lcover1:
                    nc.append(i3)
            ins = []
            while nc != []:
                id0, id1 = find_short_node(graph, pauli_map, nc, dp)
                connect_node(graph, pauli_map, pauli_map[id0], pauli_map[id1], ins)
                nc.remove(id0)
            for i3 in ins:
                if i3[0] == 'swap':
                    qc.swap(i3[1][0], i3[1][1])
            qc = dummy_synthesis([[i2]], pauli_map=pauli_map, graph=graph, qc=qc)
    return qc

def dummy_synthesis(pauli_layers, pauli_map=None, graph=None, qc=None, arch='manhattan'):
    lnq = len(pauli_layers[0][0][0]) # logical qubits
    if graph == None:
        G, C = load_graph(arch, dist_comp=True) # G is adj, C is dist
        graph = pGraph(G, C)
    if pauli_map == None:
        pauli_map = dummy_qubit_mapping(graph, lnq)
    else:
        add_pauli_map(graph, pauli_map)
    pnq = len(graph) # physical qubits
    if qc == None:
        qc = QuantumCircuit(pnq)
    for i1 in pauli_layers: # i1 is layer of blocks
        for i2 in i1: # i2 is block of pauli strings
            for i3 in i2:  # i3 is pauli string
                cns = ps2nodes(i3.ps)
                pauli_single_gates(qc, pauli_map, i3.ps, left=True)
                # for i in cns:
                #     if i3.ps[i] == 'X':
                #         qc.u(np.pi/2, 0, np.pi, pauli_map[i])
                #         # qc.h(pauli_map[i])
                #     elif i3.ps[i] == 'Y':
                #         qc.u(np.pi/2, -np.pi/2, np.pi/2, pauli_map[i])
                for i4 in range(len(cns)-1):
                    dummy_local_move(qc, graph, pauli_map, pauli_map[cns[i4]], pauli_map[cns[i4+1]])
                if len(cns) >= 1:
                    qc.rz(i3.real, pauli_map[cns[-1]])
                for i4 in range(len(cns)-1, 0, -1):
                    dummy_local_move(qc, graph, pauli_map, pauli_map[cns[i4-1]], pauli_map[cns[i4]])
                pauli_single_gates(qc, pauli_map, i3.ps, left=False)
                # for i in cns:
                #     if i3.ps[i] == 'X':
                #         # qc.h(pauli_map[i])
                #         qc.u(np.pi/2, 0, np.pi, pauli_map[i])
                #     elif i3.ps[i] == 'Y':
                #         # Y = 1/sqrt{2} [[1, i],[i, 1]]
                #         qc.u(-np.pi/2, -np.pi/2, np.pi/2, pauli_map[i])
    return qc

