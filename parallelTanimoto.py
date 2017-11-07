import csv

def open_file(chemicals):
    with open("chemicals.tsv") as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        for row in rd:
            chemicals.update({row[0]:row[2]})

def analyze_string(smiles, letters):
    for letter in smiles:
        if letter in letters.keys():
            new_value = letters.get(letter) + 1
            letters.update({letter:new_value})
        else:
            letters.update({letter:1})
    letters.update({'@':1})

def common_elements(elem_a, elem_b):
    
    
chemicals = {}
letters = {}
open_file(chemicals)
analyze_string(chemicals.get('ZINC00006923'), letters)
print(chemicals.get('ZINC00006923'))
print(letters)
