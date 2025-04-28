from random import shuffle
import random

def generate_random(qubit_num: int, single_excitation_num: int, double_excitation_num: int):
    random.seed(42)
    qubit_list_1 = list(range(qubit_num // 2))
    qubit_list_2 = list(range(qubit_num // 2, qubit_num))
    excitations = []
    for i in range(single_excitation_num):
        shuffle(qubit_list_1)
        shuffle(qubit_list_2)
        excitation = [qubit_list_1[0], qubit_list_2[0]]
        excitation.sort()
        excitations.append(excitation)
    for i in range(double_excitation_num):
        shuffle(qubit_list_1)
        shuffle(qubit_list_2)
        excitation = [qubit_list_1[0], qubit_list_1[1], qubit_list_2[0], qubit_list_2[1]]
        excitation.sort()
        excitations.append(excitation)

    shuffle(excitations)
    
    return excitations

def write_to_file(qubit_num: int, excitations: list, file_name: str):
    with open(file_name, 'a') as f:
        f.write(f'{qubit_num}\n')
        for excitation in excitations:
            for e in excitation:
                f.write(f'{e} ')
            f.write('\n')
        f.write('\n')

if __name__ == '__main__':
    for i in range(8, 40, 4):
        excitations = generate_random(i, i // 2, i // 2)
        write_to_file(i, excitations, 'benchmarks/random_tests.txt')