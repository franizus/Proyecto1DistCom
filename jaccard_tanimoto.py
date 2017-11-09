""" Parallel Jaccard/Tanimoto coefficient module """
import csv
import math
import sys
import time
from threading import Thread


def open_file():
    """Opens a tsv file and returns a list with the data."""
    chemicals_list = []
    with open("ZINC_chemicals1.tsv") as file:
        reader = csv.reader(file, delimiter="\t", quotechar='"')
        for row in reader:
            row_aux = (row[1], row[3])
            chemicals_list.append(row_aux)
    return chemicals_list


def analyze_string(chemical_compound):
    """Analyzes the string of chemical compound's elements \
    and returns the analyzed data in a dictionary."""
    analyzed_string = {}
    for character in chemical_compound:
        if character in analyzed_string:
            new_value = analyzed_string.get(character) + 1
            analyzed_string.update({character: new_value})
        else:
            analyzed_string.update({character: 1})
    if '@' in analyzed_string:
        analyzed_string.update({'@': 1})
    return analyzed_string


def number_common_elements(chemical_a, chemical_b):
    """Returns the number of common elements between chemical \
    compound a and chemical compound b."""
    copy_b = chemical_b.copy()
    number_elements = 0
    for key_a in chemical_a:
        for key_b in copy_b:
            if key_a == key_b:
                number_elements += min(chemical_a.get(key_a),
                                       chemical_b.get(key_b))
                del copy_b[key_b]
                break
    return number_elements


def number_chemical_elements(chemical_compound):
    """Returns the number of elements in a chemical compound."""
    number_elements = 0
    for key in chemical_compound:
        number_elements += chemical_compound.get(key)
    return number_elements


def jac_tan_coefficient(elements_a, elements_b, common_elements):
    """Returns the coefficient of Jaccard/Tanimoto between two chemical compounds."""
    return round((common_elements / (elements_a + elements_b - common_elements)), 2)


def get_pivots(chemicals_length, number_processors):
    """Calculates the pivots to divide the chemicals between the threads."""
    pivots_list = []
    total = (chemicals_length**2 + chemicals_length) / 2
    elements_per_processor = int(total / number_processors)
    pivots_list.append(1)
    for i in range(number_processors - 1):
        pivots_list.append(
            int(round(solve_quadratic_equation(elements_per_processor * (i + 1)))) + 1)
    pivots_list.append(chemicals_length + 1)
    return pivots_list


def solve_quadratic_equation(c_value):
    """Solves a quadratic equation in the form ax2 + bx + c = 0, given c."""
    c_value *= -2
    return abs((-1 + math.sqrt(1 - 4 * c_value)) / 2)


def fill_compared_list(chemicals_list, pivot_min, pivot_max, compared_chemicals_list):
    """Fills a list with the comparison of the chemical compounds between two pivots in the list."""
    for i in range(pivot_min, pivot_max):
        for j in range(i):
            letters_a = analyze_string(chemicals_list[i][1])
            letters_b = analyze_string(chemicals_list[j][1])
            coef = jac_tan_coefficient(number_chemical_elements(letters_a),
                                       number_chemical_elements(letters_b),
                                       number_common_elements(letters_a, letters_b))
            row_aux = (chemicals_list[i][0], chemicals_list[j][0], coef)
            compared_chemicals_list.append(row_aux)


def start_join_all(thread_array):
    """Starts and joins all the threads in the array."""
    for t in thread_array:
        t.start()
    for t in thread_array:
        t.join()


def write_file(compared_chemicals_list):
    """Writes a tsv file with the compared chemicals."""
    with open("chem_sim_total.tsv", "w") as record_file:
        record_file.write("Chem_ID_1\tChem_ID_2\tTanimoto_similarity\n")
        for row in compared_chemicals_list:
            record_file.write(row[0] + "\t" + row[1] +
                              "\t" + str(row[2]) + "\n")


if __name__ == "__main__":
    NUMBER_THREADS = int(sys.argv[1])
    CHEMICALS_LIST = open_file()
    PIVOTS_LIST = get_pivots(len(CHEMICALS_LIST) - 1, NUMBER_THREADS)

    start = time.time()
    compared_chemicals = [[], [], [], [], [], [], [], []]
    threads = []
    for index in range(NUMBER_THREADS):
        threads.append(Thread(target=fill_compared_list, args=(
            CHEMICALS_LIST, PIVOTS_LIST[index],
            PIVOTS_LIST[index + 1],
            compared_chemicals[index], )))
    start_join_all(threads)

    compared_chemicals_total = []
    for index in range(NUMBER_THREADS):
        compared_chemicals_total += compared_chemicals[index]
    
    end = time.time()
    print(end - start)
    compared_chemicals_total = sorted(
        compared_chemicals_total, key=lambda id: id[0])
    write_file(compared_chemicals_total)
