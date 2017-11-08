import csv

def open_file(chemicals):
    with open("chemicals.tsv") as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        for row in rd:
            chemicals.update({row[0]:row[2]})

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

chemicals = {}
letters_a = {}
letters_b = {}
open_file(chemicals)
analyze_string(chemicals.get('ZINC00006923'), letters_a)
analyze_string(chemicals.get('ZINC04843014'), letters_b)
print(chemicals.get('ZINC00006923'))
print(letters_a)
print(chemicals.get('ZINC04843014'))
print(letters_b)
na = number_of_elements(letters_a)
nb = number_of_elements(letters_b)
nc = common_elements(letters_a, letters_b)
coe = jac_tan_coefficient(na, nb, nc)
print(na)
print(nb)
print(nc)
print(coe)
