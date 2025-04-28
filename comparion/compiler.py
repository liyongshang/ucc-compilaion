from functions import *
from arch import *
import synthesis_SC
from parallel_bl import depth_oriented_scheduling
from .tools import print_qc
from pytket.circuit import Qubit, PauliExpBox, OpType
from qiskit import QuantumCircuit, transpile
from pytket.pauli import Pauli, QubitPauliString
from pytket import Circuit
from time import time
from pytket.passes import PauliSimp, FullPeepholeOptimise
from pytket.qasm import circuit_to_qasm_str
from pytket.transform import Transform, PauliSynthStrat, CXConfigType
from pytket.utils import gen_term_sequence_circuit, QubitPauliOperator

class Board:
    def __init__(self, graph, cs):
        self.graph = graph
        self.color = []  # 一种颜色占据的位置
        self.grid = []  # 一个位置的颜色
        self.edge = []  # 边缘位置
        # self.weight = [max(1/(1.1**l), 0.000000001) for l in layers]
        self.ct = -1
        md = 10000
        for p1 in range(len(graph.C)):
            d = 0
            for p2 in range(len(graph.C)):
                d += self.graph.C[p1][p2]
            if d < md:
                md = d
                self.ct = p1
        self.edge.append(self.ct)
        for i in range(len(graph.G)):
            self.grid.append([])
        for c in cs:
            self.color.append([])

    def move(self, c, pos):
        self.grid[pos] = c
        self.edge.remove(pos)
        for i in c:
            self.color[i].append(pos)
        for a in self.graph[pos].adj:
            if len(self.grid[a]) == 0 and a not in self.edge:
                self.edge.append(a)

    def score(self, c, p):
        td = 0
        for ci in c:
            mdci = 10000
            if len(self.color[ci]) == 0:
                mdci = 0
            for pci in self.color[ci]:
                mdci = min(mdci, self.graph.C[pci][p])
            td += mdci  # * self.weight[ci]
        td = td * 100 + self.graph.C[self.ct][p]
        return td

    def unconnected_degree(self, pss, pi):
        ud = 0
        for ps in pss:
            connected = ps[:1]
            remain = ps[1:]
            for i in range(len(ps) - 1):
                md = 10000
                nn = -1
                for q1 in connected:
                    for q2 in remain:
                        if self.graph.C[pi[q1]][pi[q2]] < md:
                            md = self.graph.C[pi[q1]][pi[q2]]
                            nn = q2
                ud += md - 1
                connected.append(nn)
                remain.remove(nn)
        return ud

class QScheduler:
    def __init__(self, bs, nq):
        pc = []
        for i in range(nq):
            pc.append([])
        i = 0
        for b in bs:
            for q in b:
                pc[q].append(i)
            i += 1
        self.pieces = []
        for p in pc:
            self.pieces.append([False] + p)
        self.placed = []

    def move(self, i):
        self.pieces[i][0] = True
        self.placed.append(i)

    def cdd_pieces(self):
        indp = []
        nq = len(self.pieces)
        for i in range(nq):
            if self.pieces[i][0] or len(self.pieces[i]) == 1:
                continue
            flag = True
            for j in range(nq):
                if j == i or self.pieces[j][0] or len(self.pieces[i]) >= len(self.pieces[j]):
                    continue
                if all(it in self.pieces[j][1:] for it in self.pieces[i][1:]):
                    flag = False
            if flag:
                indp.append(i)

        res = []
        mpriority1 = -1
        mpriority2 = -1
        for pi in indp:
            priority1 = 0
            priority2 = len(self.pieces[pi])
            for pp in self.placed:
                for c in self.pieces[pp][1:]:
                    if c in self.pieces[pi]:
                        priority1 += 1
            if priority1 > mpriority1 or (priority1 == mpriority1 and priority2 > mpriority2):
                res.clear()
                res.append(pi)
                mpriority1 = priority1
                mpriority2 = priority2
            if priority1 == mpriority1 and priority2 == mpriority2:
                res.append(pi)
        return res

class Compiler:
    def __init__(self, pauli_blocks, graph):
        self.graph = graph

        blocks = []
        for bk in pauli_blocks:
            blocks.append(compute_block_cover(bk))
        # print(blocks[:5])
        self.board = Board(self.graph, blocks)
        self.scheduler = QScheduler(blocks, len(pauli_blocks[0][0]))

        self.lnq = len(pauli_blocks[0][0])
        self.pauli_layers = [[b] for b in pauli_blocks]
        self.pauli_layers_ph = depth_oriented_scheduling(pauli_blocks, length=self.lnq // 2, maxiter=30)
        self.op_list = pauli_blocks

    def set_phycir_path(self, path):
        self.phycir_path = path
    
    def ph_compile(self, opt=0):
        self.my_pi = dummy_mapping(self.lnq)
        if opt==1:
            self.initial_mapping()
        qc0 = QuantumCircuit(self.lnq)
        # qc0.x([self.my_pi[0]])
        qc, inner, outer = synthesis_SC.block_opt_SC(self.pauli_layers, graph=self.graph, pauli_map=self.my_pi)
        # f = open(self.phycir_path + '_origin.txt', mode='w+')
        # f.write(qc.qasm())
        # f.close()
        # ncx, nsg, qc = print_qc(qc, opt_level=3)
        # f = open(self.phycir_path + '.txt', mode='w+')
        # f.write(qc.qasm())
        # f.close()
        # print('ph result: inner: ', inner, ': outer: ', outer, '\n')
        # gate_names = set()
        # for instr, qargs, cargs in qc.data:
        #     gate_names.add(instr.name)
        # qc = transpile(qc, basis_gates=list(gate_names), optimization_level=1)
        return qc # [0, 0, inner + outer, outer, qc.depth()]

    def initial_mapping(self):
        self.my_pi = {}
        while True:
            cdd_p = self.scheduler.cdd_pieces()
            if len(cdd_p) == 0:
                break
            c = -1
            pos = -1
            mscore = 1000000000
            for posi in self.board.edge:
                for ci in cdd_p:
                    score = self.board.score(self.scheduler.pieces[ci][1:], posi)
                    if score < mscore:
                        mscore = score
                        c = ci
                        pos = posi
            self.scheduler.move(c)
            self.board.move(self.scheduler.pieces[c][1:], pos)
            self.my_pi[c] = pos
            print(c, ' ', end="")
        print('\n')
        return self.my_pi
    
    def go_compile(self, opt=0):
        self.initial_mapping()
        if opt == 1:
            self.my_pi = dummy_mapping(self.lnq)
        qc0 = QuantumCircuit(self.lnq)
        qc0.x([self.my_pi[0], self.my_pi[1]])
        qc, inner, outer = synthesis_SC.block_opt_SC(self.pauli_layers_ph, graph=self.graph, pauli_map=self.my_pi, synthesis_opt=True)
        # f = open(self.phycir_path + '_opt' + str(opt) + '_origin.txt', mode='w+')
        # f.write(qc.qasm())
        # f.close()
        # ncx, nsg, qc = print_qc(qc, opt_level=3)
        # f = open(self.phycir_path + '_opt' + str(opt) + '.txt', mode='w+')
        # f.write(qc.qasm())
        # f.close()
        # print('my result: inner: ', inner, ', outer: ', outer, ', depth: ', qc.depth())
        return qc0.compose(qc, list(range(self.lnq))) # [0, 0, inner + outer, outer, qc.depth()]
    
    def start(self, cp, opt=0):
        if cp == 'ph':
            return self.ph_compile(opt=opt)
        if cp == 'go':
            return self.go_compile(opt=opt)
        if cp == 'tk':
            return self.tk_ucc_compile()
        
    def tk_compile(self):
        def to_pauli_list(ps):
            r = []
            for i in ps:
                if i == 'I':
                    r.append(Pauli.I)
                elif i == 'X':
                    r.append(Pauli.X)
                elif i == 'Y':
                    r.append(Pauli.Y)
                elif i == 'Z':
                    r.append(Pauli.Z)
            return r
        coupling = []
        for i in range(len(self.graph.G)):
            for j in range(len(self.graph.G)):
                if self.graph.G[i][j] != 0:
                    coupling.append([i,j])
        oplist = {}
        n = len(self.op_list[0][0])
        q = [Qubit(i) for i in range(n)]
        for i in self.op_list:
            for j in i:
                op = QubitPauliString(q, to_pauli_list(j.ps))
                oplist[op] = 1/3.14
        def add_excitation(circ, term_dict, param=1.0):
            for term, coeff in term_dict.items():
                qubits, paulis = zip(*term.map.items())
                pbox = PauliExpBox(paulis, coeff * param)
                circ.add_pauliexpbox(pbox, qubits)
        ansatz = Circuit(n)
        t0 = time()
        add_excitation(ansatz, oplist)
        PauliSimp().apply(ansatz)
        # print(f"TK Pauli Simp: {time()-t0}")
        # t0 = time()
        # FullPeepholeOptimise().apply(ansatz)
        # print(f"TK O2: {time()-t0}")
        # t0 = time()
        qstr = circuit_to_qasm_str(ansatz)
        qc = QuantumCircuit.from_qasm_str(qstr)
        # t0 = time()
        qc = transpile(qc, basis_gates=['cx', 'swap', 'u3'], coupling_map=coupling)
        # print("Qiskit L3:", time()-t0)
        return qc
    
    def tk_ucc_compile(self):
        qubit_list = [Qubit(i) for i in range(self.lnq)]
        qps_dict = {}
        for b in self.op_list:
            for pauli_str in b:
                temp_string = []
                for pauli in pauli_str.ps:
                    if pauli == 'X':
                        temp_string.append(Pauli.X)
                    elif pauli == 'Y':
                        temp_string.append(Pauli.Y)
                    elif pauli == 'Z':
                        temp_string.append(Pauli.Z)
                    else:
                        temp_string.append(Pauli.I)
                # qps_list.append(QubitPauliString(qubit_list, temp_string))
                qps_dict[QubitPauliString(qubit_list, temp_string)] = pauli_str.coeff

        operator = QubitPauliOperator(qps_dict)
        init_circ = Circuit(len(operator.all_qubits))
        set_synth_circuit = gen_term_sequence_circuit(
            operator, init_circ
        )
        Transform.UCCSynthesis(
        PauliSynthStrat.Sets, CXConfigType.Tree
        ).apply(set_synth_circuit)
        set_synth_cx_count = set_synth_circuit.n_gates_of_type(OpType.CX)
        set_synth_cx_depth = set_synth_circuit.depth_by_type(OpType.CX)

        # print(f"\tAfter Set synth: CX count: {set_synth_cx_count}, depth: {set_synth_cx_depth}")

        qstr = circuit_to_qasm_str(set_synth_circuit)
        qstr = qstr.replace('*I', '')
        # f = open('./data/debug.qasm', 'w+')
        # f.write(qstr)
        # f.close()
        qc = QuantumCircuit.from_qasm_str(qstr)
        return qc