from base import *
import pickle
from qiskit import QuantumCircuit
import numpy
import steiner_forest

class para:
    def __init__(self, name, coeff) -> None:
        self.name = name
        self.coeff = coeff

def interaction(q, Q, block_size):
    res = 0
    for a in Q:
        for i in range(len(q)):
            res += q[i] * a[i] / (block_size[i] ** 2)
    return res

def rearrange(Q):
    block_size = [0 for i in range(len(Q[0]))]
    for i in range(len(block_size)):
        for j in range(len(Q)):
            block_size[i] += Q[j][i]
    res = []
    permutation = []
    rows = list(range(len(Q)))
    iMaxWeightRow = 0
    for i in rows:
        if sum(Q[i]) > sum(Q[iMaxWeightRow]):
            iMaxWeightRow = i
    res.append(Q[iMaxWeightRow])
    rows.remove(iMaxWeightRow)
    permutation.append(iMaxWeightRow)
    while (len(rows)):
        iInteractMax = rows[0]
        for i in rows:
            if interaction(Q[i], res, block_size) > interaction(Q[iInteractMax], res, block_size):
                iInteractMax = i
        res.append(Q[iInteractMax])
        rows.remove(iInteractMax)
        permutation.append(iInteractMax)
    return res, permutation

def simulation(s, Q, order): # pi: map<int,int>; Q: vector<vector<int>>; order: vector<int>; C: vector<vector<int>>
    m = s.m.copy()
    return steiner_forest.simulation(m)
    nMappedQubits = len(m)
    for i in order[nMappedQubits:]:
        q = Q[i]
        md = maxsize
        bestVertex = -1
        for vertex in range(nHardwareQubit):
            d = 0
            if vertex in m.values():
                continue
            for k in m.keys():
                _q = Q[k]
                d += sum([a * b for a, b in zip(q, _q)]) * C[vertex][m[k]]
            if d < md:
                md = d
                bestVertex = vertex
        m[i] = bestVertex
    return m

def showMappingOnGrid(m, a, b):
    print('------------------------')
    matrix = [['+' for i in range(b)] for j in range(a)]
    for k, v in m.items():
        i = v // b
        j = v % b
        matrix[i][j] = k
    max_width = max(len(str(element)) for row in matrix for element in row) + 2

    # Print the matrix with alignment
    for row in matrix:
        print(" ".join(f"{element:>{max_width}}" for element in row))

def count_swap(gate_list):
    count = 0
    for gate in gate_list:
        if gate[0] == "SWAP":
            count += 1
    return count


def gate_to_str(gate):
    result = ""
    for g in gate:
        if isinstance(g, para):
            result += f"{g.name}*{g.coeff} "
        else:
            result += f"{g} "
    return result

from qiskit.circuit import Parameter
def trans_qiskit(gate_list: list):
    parameter = Parameter('phi')
    qc = QuantumCircuit(nHardwareQubit)
    for gate in gate_list:
        gate_name = gate[0]
        if gate_name == "SWAP":
            qc.swap(gate[1][0], gate[1][1])
        elif gate_name == "CX":
            qc.cx(gate[1][0], gate[1][1])
        elif gate_name == "CZ":
            qc.cz(gate[1][0], gate[1][1])
        elif gate_name == "RY":
            param = gate[2]
            if isinstance(param, para):
                qc.ry(parameter * param.coeff * numpy.pi, gate[1][0])
            else:
                qc.ry(gate[2], gate[1][0])
        elif gate_name == "RZ":
            param = gate[2]
            if isinstance(param, para):
                qc.rz(param, gate[1][0])
            else:
                qc.rz(gate[2], gate[1][0])
        elif gate_name == "H":
            qc.h(gate[1][0])
        elif gate_name == "X":
            qc.x(gate[1][0])
        else:
            raise ValueError("Invalid gate name")
    return qc

def count_excitation(benchmark):
    with open(f"./benchmarks/{benchmark}.pickle", "rb") as f:
        blocks = pickle.load(f)

    single_excitation, double_excitation = 0, 0
    for block in blocks:
        ps = block[0][0]
        cnt = 0
        for item in ps:
            if item in ["X", "Y"]:
                cnt += 1
        if cnt == 2:
            single_excitation += 1
        elif cnt == 4:
            double_excitation += 1
        else:
            raise ValueError("Invalid excitation number")

    with open(f"./excitations/overall_excitation.txt", "a") as f:
        f.write(f"{benchmark}: \n")
        f.write(f"single excitation: {single_excitation}\n")
        f.write(f"double excitation: {double_excitation}\n\n")
    with open(f"./excitations/{benchmark}_excitation.txt", "w") as f:
        f.write(f"{benchmark}: \n")
        f.write(f"single excitation: {single_excitation}\n")
        f.write(f"double excitation: {double_excitation}\n")
        f.write(f"blocks: \n")
        for block in blocks:
            f.write(f"{block}\n")


def reverse_ps(benchmark: str):
    with open(f"./benchmarks/{benchmark}.pickle", "rb") as f:
        blocks = pickle.load(f)

    for block in blocks:
        for b in block:
            b[0] = b[0][::-1]

    pickle.dump(blocks, open(f"./benchmarks_reverse/{benchmark}.pickle", "wb"))