""" Parallel Jaccard/Tanimoto coefficient module """
import csv
import math
import time
import multiprocessing
from threading import Thread


def open_file():
    """Opens a tsv file and returns a list with the data."""
    chemicals_list = []
    with open("ZINC_chemicals.tsv") as file:
        reader = csv.reader(file, delimiter="\t", quotechar='"')
        for row in reader:
            row_aux = (row[1], row[3])
            chemicals_list.append(row_aux)
    chemicals_list = sorted(chemicals_list, key=lambda id: id[0])
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
    tolal_elements = (chemicals_length**2) / 2
    pivots_list.append(0)
    for i in range(1, number_processors):
        pivots_list.append(
            int(round(tolal_elements - (math.sqrt(i / number_processors) * chemicals_length))))
    pivots_list.append(chemicals_length - 1)
    return pivots_list


def fill_compared_list(chemicals_list, pivot_min, pivot_max, compared_chemicals_list):
    """Fills a list with the comparison of the chemical compounds between two pivots in the list."""
    for i in range(pivot_min, pivot_max):
        for j in range(i + 1, len(chemicals_list)):
            letters_a = analyze_string(chemicals_list[i][1])
            letters_b = analyze_string(chemicals_list[j][1])
            coef = jac_tan_coefficient(number_chemical_elements(letters_a),
                                       number_chemical_elements(letters_b),
                                       number_common_elements(letters_a, letters_b))
            row = (chemicals_list[i][0], chemicals_list[j][0], coef)
            compared_chemicals_list.append(row)


def start_join_all(thread_array):
    """Starts and joins all the threads in the array."""
    for thread in thread_array:
        thread.start()
    for thread in thread_array:
        thread.join()


def print_to_console(number_threads, compared_chemicals_list):
    """Prints in console the compared chemicals."""
    for i in range(number_threads):
        for j in compared_chemicals_list[i]:
            print(j)


if __name__ == "__main__":
    NUMBER_THREADS = multiprocessing.cpu_count()
    CHEMICALS_LIST = open_file()
    PIVOTS_LIST = get_pivots(len(CHEMICALS_LIST), NUMBER_THREADS)

    START_TIME = time.time()
    compared_chemicals = []
    for index in range(NUMBER_THREADS):
        compared_chemicals.append([])
    threads = []
    for index in range(NUMBER_THREADS):
        threads.append(Thread(target=fill_compared_list, args=(
            CHEMICALS_LIST, PIVOTS_LIST[index],
            PIVOTS_LIST[index + 1],
            compared_chemicals[index], )))
    start_join_all(threads)

    END_TIME = time.time()
    print_to_console(NUMBER_THREADS, compared_chemicals)
    print(END_TIME - START_TIME)
