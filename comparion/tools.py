def calc_qc(qc):
    from qiskit import transpile
    # qc = transpile(qc, basis_gates=['cx', 'u3'], optimization_level=0)
    c = qc.count_ops()
    t0 = sum(c.values())
    # print(c)
    t2 = 0
    t1 = 0
    if 'cx' in c:
        t1 = c['cx']
    if 'swap' in c:
        t2 = c['swap']
    return t1+t2*3, t0-t1-t2, qc.depth()
    
def count_sched(sched):
    c = 0
    s = 0
    nq = len(sched[0][0][0]) # sched : pauli layers
    ns = 0
    for i in sched:
        for k in i:
            for j in k:
                c += 2*max(nq - 1- j.count('I'), 0)
                s += 2*(nq - j.count('I') - j.count('Z')) + 1 # accurate
                ns += 1
    return nq, ns, c, s
def count_oplist(parr):
    sched = [[i] for i in parr]
    return count_sched(sched)

import sys
from qiskit import transpile
def print_qc(qc, f=sys.stdout, opt_level=0):
    qc = transpile(qc, basis_gates=['cx', 'u3'], optimization_level=opt_level)
    c = qc.count_ops()
    t0 = sum(c.values())
    if 'cx' in c:
        t1 = c['cx']
    else:
        t1 = 0
    if f != None:
        print('CNOT: ' + str(t1) + ", Single: " + str(t0-t1) + ', Total: ' + str(t0) + ', Depth:', qc.depth(), file=f, flush=True)
    return t1, t0-t1, qc
    
import os    
def set_cwd():
    def get_script_path():
        return os.path.dirname(os.path.realpath(sys.argv[0]))
    os.chdir(get_script_path())

from mypauli import pauliString
from copy import deepcopy
import numpy as np
def fermi2ps(fermi_ops, lnq):
    xy1 = ['XY', 'YX']
    sig1 = [-np.pi * 0.5, np.pi * 0.5]
    xy2 = ['XYXX', 'YXXX', 'YYYX', 'YYXY', 'XXYX', 'XXXY', 'YXYY', 'XYYY']
    sig2 = [-np.pi * 0.125] * 4 + [np.pi * 0.125] * 4
    pauli_blocks = []
    for op in fermi_ops:
        b = []
        op.sort()
        if len(op) == 2:
            for i in range(2):
                ps = ['I'] * lnq
                ps[op[0]] = xy1[i][0]
                ps[op[1]] = xy1[i][1]
                for j in range(op[0] + 1, op[1]):
                    ps[j] = 'Z'
                b.append(pauliString(''.join(ps), sig1[i]))
            pauli_blocks.append(deepcopy(b))
        if len (op) == 4:
            for i in range(8):
                ps = ['I'] * lnq
                for j in range(4):
                    ps[op[j]] = xy2[i][j]
                for j in list(range(op[0] + 1, op[1])) + list(range(op[2] + 1, op[3])):
                    ps[j] = 'Z'
                b.append(pauliString(''.join(ps), sig2[i]))
            pauli_blocks.append(deepcopy(b))
    return pauli_blocks
                