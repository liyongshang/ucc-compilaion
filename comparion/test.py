from compiler import Compiler
from mypauli import *
import pickle


def program_prep(origin_blocks):
    bn = pn = 0
    blocks = []
    for bk in origin_blocks:
        blocks.append([pauliString(ps[0], ps[1]) for ps in bk])
    bn = len(origin_blocks)
    return blocks
def compile(benchmark):
    with open(f'../benchmarks/{benchmark}.pickle', 'rb') as f:
        op_list = pickle.load(f)
    blocks = program_prep(op_list)
    compiler = Compiler(blocks, 'grid0505')
    qc = compiler.start('ucc')
    # print(qc.count_ops())
    print(qc.count_ops())
    print(qc.depth())
    # qc.draw('mpl',filename='./result/LiH.png')
    f = open(f'./result/{benchmark}.qasm', 'w+')
    f.write(qc.qasm())
    f.close()

benchmarks = ["LiH", "HF", "BeH2", "H2O", "NH3", "CH4", "NaH", "N2"]
# benchmarks = ['LiH']
for benchmark in benchmarks:
    compile(benchmark)

