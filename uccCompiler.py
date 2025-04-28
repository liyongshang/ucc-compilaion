from MCTS import MCTS
from base import *
from tools import *
from evaluation import addDoubleFermionic
from time import time
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), "comparion"))
from evaluation import evalPaulis
from tools import simulation
from comparion.compiler import Compiler
from comparion.tools import fermi2ps


def startComp(excitations, cp: str, lnq):
    pauli_blocks = fermi2ps(excitations, lnq)
    comp = Compiler(pauli_blocks, Graph)
    qc = comp.start(cp)
    return qc

def femi2Q(femi_ops, nQubits):
    B = []
    for op in femi_ops:
        op.sort()
        B.append(op.copy())
        if len(op) == 2:
            B.append(op + list(range(op[0]+1, op[1])))
        if len(op) == 4:
            B.append(op + list(range(op[0]+1, op[1])) + list(range(op[2]+1, op[3])))
    Q = [[0 for _ in B] for i in range(nQubits)]
    for j in range(len(B)):
        for i in B[j]:
            Q[i][j] = 1
    return Q

def start(excitations, cp, lnq, retPi = False):
    Q = femi2Q(excitations, lnq)
    myMCTS = MCTS(Q, evalPaulis, simulation, excitations, Graph.find_center())
    # print('challenge times: ', mChallengeTimes)
    myMCTS.setChallengeTimes(mChallengeTimes)
    pi = {i: i for i in range(lnq)}

    if cp == 'heu':
        pi = myMCTS.heuristicMapping()

    if cp == 'mc':
        t0 = time()
        myMCTS.start()
        cp_time = time() - t0
        print('MCTS compile time: {}'.format(cp_time))
        pi = myMCTS.tree.s.m

    conAware = True
    if cp == 'ss':
        conAware = False
    gate_list_all = []
    for op in excitations:
        op.sort()
        gate_list = []
        addDoubleFermionic(Graph, op, gate_list, conAware, pi)
        gate_list_all.append(gate_list)

    gates = []
    for gate_list in gate_list_all:
        gates.extend(gate_list)
    qc = trans_qiskit(gates)
    if retPi:
        return qc, pi
    else:
        return qc

from qiskit import transpile
def compile(excitations, cp, lnq, retPi = False):
    if cp == 'tk+SABRE':
        qc = startComp(excitations, 'tk', lnq)
        coupling = []
        for i in range(len(Graph.G)):
            for j in range(len(Graph.G)):
                if Graph.G[i][j] != 0:
                    coupling.append([i,j])
        gate_names = set()
        for instr, qargs, cargs in qc.data:
            gate_names.add(instr.name)
        qc = transpile(qc, basis_gates=list(gate_names) + ['swap'], coupling_map=coupling,\
                        layout_method='sabre', routing_method='sabre', optimization_level=0)
        return qc
    if cp in ['tk', 'ph']:
        return startComp(excitations, cp, lnq)
    else:
        return start(excitations, cp, lnq, retPi)