from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit.library.standard_gates import SwapGate

from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit.library.standard_gates import SwapGate

from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit.library.standard_gates import SwapGate

def get_ancilla_qubits_number(physical_circuit: QuantumCircuit, initial_mapping: dict) -> int:
    # 第一步：收集所有被使用的物理比特
    used_physical_qubits = set()
    for instr, qargs, _ in physical_circuit.data:
        for qubit in qargs:
            used_physical_qubits.add(physical_circuit.find_bit(qubit).index)

    # 第二步：确定逻辑比特数和辅助比特数
    logical_qubits = list(initial_mapping.keys())
    mapped_physical_qubits = set(initial_mapping.values())
    ancilla_physical_qubits = used_physical_qubits - mapped_physical_qubits
    return len(ancilla_physical_qubits)

def recover_logical_circuit(physical_circuit: QuantumCircuit, initial_mapping: dict, ancilla_qubits_number) -> QuantumCircuit:
    # 第一步：收集所有被使用的物理比特
    used_physical_qubits = set()
    for instr, qargs, _ in physical_circuit.data:
        for qubit in qargs:
            used_physical_qubits.add(physical_circuit.find_bit(qubit).index)

    # 第二步：确定逻辑比特数和辅助比特数
    logical_qubits = list(initial_mapping.keys())
    mapped_physical_qubits = set(initial_mapping.values())
    ancilla_physical_qubits = used_physical_qubits - mapped_physical_qubits

    # 给逻辑比特和辅助比特分别编号
    logical_qreg = QuantumRegister(len(logical_qubits), name='q')
    ancilla_qreg = QuantumRegister(ancilla_qubits_number, name='anc')

    logical_circuit = QuantumCircuit(logical_qreg, ancilla_qreg)

    # 建立物理比特 -> (量子寄存器, 索引) 的映射
    physical_to_logical = {}  # 物理比特 → (寄存器, 索引)
    logical_to_physical = {}

    # 逻辑比特
    for logical_idx, physical_idx in initial_mapping.items():
        physical_to_logical[physical_idx] = (logical_qreg, logical_idx)
        logical_to_physical[logical_idx] = physical_idx

    # 辅助比特
    ancilla_list = sorted(ancilla_physical_qubits)  # 排序保证编号稳定
    for ancilla_idx, physical_idx in enumerate(ancilla_list):
        physical_to_logical[physical_idx] = (ancilla_qreg, ancilla_idx)

    # 第三步：遍历指令，重建电路
    for instr, qargs, cargs in physical_circuit.data:
        physical_indices = [physical_circuit.find_bit(qubit).index for qubit in qargs]

        if isinstance(instr, SwapGate):
            q0, q1 = physical_indices
            l0 = physical_to_logical.get(q0, None)
            l1 = physical_to_logical.get(q1, None)

            if l0 is not None and l1 is not None:
                # 交换
                physical_to_logical[q0], physical_to_logical[q1] = l1, l0
                # 如果是逻辑比特，还要更新 logical_to_physical
                if l0[0] == logical_qreg:
                    logical_to_physical[l0[1]] = q1
                if l1[0] == logical_qreg:
                    logical_to_physical[l1[1]] = q0
            elif l0 is not None and l1 is None:
                # l0移动到q1
                physical_to_logical[q1] = l0
                del physical_to_logical[q0]
                if l0[0] == logical_qreg:
                    logical_to_physical[l0[1]] = q1
            elif l0 is None and l1 is not None:
                # l1移动到q0
                physical_to_logical[q0] = l1
                del physical_to_logical[q1]
                if l1[0] == logical_qreg:
                    logical_to_physical[l1[1]] = q0
            else:
                # 两边都是空闲比特，不处理
                pass
        else:
            # 普通门：映射物理比特
            logical_qargs = []
            for p in physical_indices:
                if p not in physical_to_logical:
                    raise ValueError(f"Physical qubit {p} has no logical or ancilla mapping.")
                reg, idx = physical_to_logical[p]
                logical_qargs.append(reg[idx])
            logical_circuit.append(instr, logical_qargs, cargs)

    return logical_circuit

import numpy as np
from qulacs import QuantumCircuit as QulacsCircuit, QuantumState
from qiskit.circuit.library import U3Gate, CXGate

def qiskit_to_qulacs(qc: QuantumCircuit) -> QulacsCircuit:
    """将 Qiskit QuantumCircuit 转换为 Qulacs QuantumCircuit。"""
    n = qc.num_qubits
    qulacs_circuit = QulacsCircuit(n)
    
    for instr, qargs, cargs in qc.data:
        if isinstance(instr, U3Gate):
            theta, phi, lam = instr.params
            qulacs_circuit.add_U3_gate(qargs[0].index, theta, phi, lam)
        elif isinstance(instr, CXGate):
            qulacs_circuit.add_CNOT_gate(qargs[0].index, qargs[1].index)
        else:
            raise NotImplementedError(f"Unsupported gate: {instr.name}")
    
    return qulacs_circuit

def are_circuits_equivalent(qc1: QuantumCircuit, qc2: QuantumCircuit, tol=1e-6) -> bool:
    """判断两个 Qiskit 电路是否等价（忽略全局相位）。"""
    if qc1.num_qubits != qc2.num_qubits:
        raise ValueError("Circuits must have the same number of qubits!")

    n = qc1.num_qubits

    # 初始化 |0...0> 初态
    state1 = QuantumState(n)
    state2 = QuantumState(n)

    # 转成 qulacs 电路
    circuit1 = qiskit_to_qulacs(qc1)
    circuit2 = qiskit_to_qulacs(qc2)

    # 应用电路
    circuit1.update_quantum_state(state1)
    circuit2.update_quantum_state(state2)

    # 获取最终 statevector
    vec1 = state1.get_vector()
    vec2 = state2.get_vector()

    # 计算 Fidelity
    fidelity = abs(np.vdot(vec1, vec2)) ** 2

    return abs(1 - fidelity) < tol

from uccCompiler import compile
import pickle
from qiskit import transpile
from generate_random import generate_random
from base import *
from time import time

def verify():
    molecules = ["LiH", "HF", "BeH2", "H2O", "NH3", "CH4", "NaH", "N2"]
    testPrograms = pickle.load(open('./benchmarks/uccs.pickle', 'rb'))
    qubitNumbers = [12 , 12 , 14 , 14 , 16 , 18 , 20 , 20]
    cps = ['ph', 'mc']
    testId = 0
    for excitations, lnq in zip(testPrograms[:8], qubitNumbers):
        # break
        print(molecules[testId])
        # excitations = excitations[:7]
        qc1 = compile(excitations, 'tk', lnq)
        # qc1 = transpile(qc1, basis_gates=['u3', 'cx'], optimization_level=3)
        qc1 = qc1.assign_parameters([-1] * len(qc1.parameters))
        pi1 = {i: i for i in range(lnq)}
        print('Paulihedral depth: ', qc1.depth())

        qc2, pi2 = compile(excitations, 'mc', lnq, True)
        # qc2 = transpile(qc2, basis_gates=['u3', 'cx'], optimization_level=3)
        qc2 = qc2.assign_parameters([1] * len(qc2.parameters))
        print('Our MCTS depth: ', qc2.depth())

        anq = max(get_ancilla_qubits_number(qc1, pi1), get_ancilla_qubits_number(qc2, pi2))
        # print(anq)
        qc1 = recover_logical_circuit(qc1, pi1, anq)
        qc2 = recover_logical_circuit(qc2, pi2, anq)
        qc1 = transpile(qc1, basis_gates=['u3', 'cx'], optimization_level=3)
        qc2 = transpile(qc2, basis_gates=['u3', 'cx'], optimization_level=3)
        t0 = time()
        eq = are_circuits_equivalent(qc1, qc2)
        print('verification time: ', time() - t0, ' s')
        print('is equavalent? ', eq)
        print('-------------')
        testId += 1

    for lnq in list(range(8, 22, 2)):
        print(f'random: {lnq}')
        excitations = generate_random(lnq, lnq // 2, lnq //2)
        excitations = excitations[:2]
        print(excitations)
        qc1 = compile(excitations, 'tk', lnq)
        # qc1 = transpile(qc1, basis_gates=['u3', 'cx'], optimization_level=3)
        qc1 = qc1.assign_parameters([-1] * len(qc1.parameters))
        pi1 = {i: i for i in range(lnq)}
        print('Paulihedral depth: ', qc1.depth())

        qc2, pi2 = compile(excitations, 'mc', lnq, True)
        # qc2 = transpile(qc2, basis_gates=['u3', 'cx'], optimization_level=3)
        qc2 = qc2.assign_parameters([1] * len(qc2.parameters))
        print('Our MCTS depth: ', qc2.depth())

        anq = max(get_ancilla_qubits_number(qc1, pi1), get_ancilla_qubits_number(qc2, pi2))
        # print(anq)
        qc1 = recover_logical_circuit(qc1, pi1, anq)
        qc2 = recover_logical_circuit(qc2, pi2, anq)
        qc1 = transpile(qc1, basis_gates=['u3', 'cx'], optimization_level=3)
        qc2 = transpile(qc2, basis_gates=['u3', 'cx'], optimization_level=3)
        t0 = time()
        eq = are_circuits_equivalent(qc1, qc2)
        print('verification time: ', time() - t0, ' s')
        print('is equavalent? ', eq)
        print('-------------')
        testId += 1

def twoCircuit():
    qc1 = QuantumCircuit(2)
    qc1.cx(0, 1)

    qc2 = QuantumCircuit(2)
    qc2.h(0) 
    qc2.h(1)
    qc2.cx(1, 0)
    qc2.h(0) 
    qc2.h(1)

    eq = are_circuits_equivalent(qc1, qc2)
    print('is equavalent? ', eq)

# ========== 使用示例 ==========

if __name__ == "__main__":
    # twoCircuit()
    verify()