import pickle
from qiskit import transpile
from generate_random import generate_random
from base import *
from uccCompiler import compile

def table1():
    molecules = ["LiH", "HF", "BeH2", "H2O", "NH3", "CH4", "NaH", "N2"]
    testPrograms = pickle.load(open('./benchmarks/uccs.pickle', 'rb'))
    qubitNumbers = [12 , 12 , 14 , 14 , 16 , 18 , 20 , 20]
    cps = ['ph', 'tk+SABRE', 'heu', 'mc']
    testId = 0
    global mapping
    for excitations, lnq in zip(testPrograms[:8], qubitNumbers):
        print(molecules[testId])
        for cp in cps:
            qc = compile(excitations, cp, lnq)
            qc = transpile(qc, basis_gates=['u3', 'cx'], optimization_level=3)
            cnot_count = 0 if 'cx' not in qc.count_ops().keys() else qc.count_ops()['cx']
            print(f'{cp}: ', cnot_count, ' ', qc.depth())
        print('-------------')
        testId += 1

    for lnq in list(range(8, 49, 2)):
        excitations = generate_random(lnq, lnq // 2, lnq //2)
        print('random', lnq)
        for cp in cps:
            qc = compile(excitations, cp, lnq)
            qc = transpile(qc, basis_gates=['u3', 'cx'], optimization_level=3)
            cnot_count = 0 if 'cx' not in qc.count_ops().keys() else qc.count_ops()['cx']
            print(f'{cp}: ', cnot_count, ' ', qc.depth())
        print('-------------')
        testId += 1

if __name__ == "__main__":
    table1()