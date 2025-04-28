import numpy as np
from pyscf import gto, scf, cc

def get_truncated_excitations(molecule, basis='sto-3g', threshold=1e-6):
    """
    计算 UCCSD 的单/双激发算符，并根据振幅筛选高贡献度的算符（不考虑自旋）。
    :param molecule: 分子结构字符串，例如 'H 0 0 0; F 0 0 0.917'
    :param basis: 计算基组，默认为 'sto-3g'
    :param threshold: 截断阈值，默认 1e-2
    :return: (单激发列表, 双激发列表)
    """
    # 构建分子
    mol = gto.M(atom=molecule, basis=basis, spin=0, charge=0)
    mf = scf.RHF(mol).run()
    mycc = cc.CCSD(mf).run()
    
    n_mo = mol.nao  # 轨道总数
    n_occ = mol.nelec[0]  # 占据轨道数
    print(n_mo, ' ', n_occ)
    
    t1 = mycc.t1  # 单激发振幅 (n_occ, n_vir)
    print(t1)
    exit()
    t2 = mycc.t2  # 双激发振幅 (n_occ, n_occ, n_vir, n_vir)
    print(t2)
    exit()
    
    single_excitations = []
    double_excitations = []
    
    # 筛选单激发
    for i in range(n_occ):
        for a in range(n_occ, n_mo):
            if abs(t1[i, a - n_occ]) > threshold:
                single_excitations.append((i, a))
    
    # 筛选双激发
    for i in range(n_occ):
        for j in range(i+1, n_occ):
            for a in range(n_occ, n_mo):
                for b in range(a+1, n_mo):
                    if abs(t2[i, j, a - n_occ, b - n_occ]) > threshold:
                        double_excitations.append((i, j, a, b))
    
    return single_excitations, double_excitations

if __name__ == "__main__":
    molecule = 'O 0 0 0; H 0 0 0.96; H 0 0 1.92'  # 例如 HF 分子
    single_exc, double_exc = get_truncated_excitations(molecule)
    
    print("筛选后的单激发算符:", single_exc)
    print("筛选后的双激发算符:", double_exc)
