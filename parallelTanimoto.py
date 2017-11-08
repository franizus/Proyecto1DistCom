import csv
import math

def open_file(chemicals):
    with open("chemicals.tsv") as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        for row in rd:
            row_aux = (row[0], row[2])
            chemicals.append(row_aux)

def analyze_string(smiles, letters):
    for letter in smiles:
        if letter in letters:
            new_value = letters.get(letter) + 1
            letters.update({letter:new_value})
        else:
            letters.update({letter:1})
    if '@' in letters:
        letters.update({'@':1})

def common_elements(elem_a, elem_b):
    copy_b = elem_b.copy()
    nc = 0
    for key_a in elem_a:
        for key_b in copy_b:
            if key_a == key_b:
                nc += min(elem_a.get(key_a), elem_b.get(key_b))
                del copy_b[key_b]
                break
    return nc

def number_of_elements(elem):
    number = 0
    for key in elem:
        number += elem.get(key)
    return number

def jac_tan_coefficient(na, nb, nc):
    return round((nc / (na + nb - nc)), 2)

def get_pivots(n, p, pivots):
    total = (n**2 + n) / 2
    elem_per_div = int(total / p)
    for i in range(p - 1):
        pivots.append(int(round(solve_quadratic(elem_per_div * (i + 1)))) + 1)

def solve_quadratic(n):
    n *= -2
    return abs((-1 + math.sqrt(1 - 4*n)) / 2)

def cal_coef(chemicals, pivot_min, pivot_max, chem_sim):
    for i in range(pivot_min, pivot_max):
        for j in range(i):
            letters_a = {}
            letters_b = {}
            analyze_string(chemicals[i][1], letters_a)
            analyze_string(chemicals[j][1], letters_b)
            coef = jac_tan_coefficient(number_of_elements(letters_a), 
                    number_of_elements(letters_b), 
                    common_elements(letters_a, letters_b))
            row_aux = (chemicals[i][0], chemicals[j][0], coef)
            chem_sim.append(row_aux)


chemicals = []
letters_a = {}
letters_b = {}
open_file(chemicals)
#chemicals = sorted(chemicals, key=lambda id: id[0])
analyze_string(chemicals[1][1], letters_a)
analyze_string(chemicals[4][1], letters_b)
print(chemicals[1][1])
print(letters_a)
print(chemicals[4][1])
print(letters_b)
na = number_of_elements(letters_a)
nb = number_of_elements(letters_b)
nc = common_elements(letters_a, letters_b)
coe = jac_tan_coefficient(na, nb, nc)
print(na)
print(nb)
print(nc)
print(coe)

pivots = []
get_pivots(len(chemicals) - 1, 4, pivots)
print(pivots)
chem_sim = []
cal_coef(chemicals, pivots[0], pivots[1], chem_sim)
chem_sim = sorted(chem_sim, key=lambda id: id[0])
print(chem_sim)
