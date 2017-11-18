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


def get_number_common_elements(chemical_a, chemical_b):
    """Returns the number of common elements between chemical \
    compound a and chemical compound b."""
    number_elements = 0
    for key_a in chemical_a:
        if key_a in chemical_b:
            number_elements += min(chemical_a.get(key_a), chemical_b.get(key_a))
    return number_elements


def get_number_chemical_elements(chemical_compound):
    """Returns the number of elements in a chemical compound."""
    number_elements = 0
    for key in chemical_compound:
        number_elements += chemical_compound.get(key)
    return number_elements


def get_jac_tan_coefficient(chemical_a, chemical_b):
    """Returns the coefficient of Jaccard/Tanimoto between two chemical compounds."""
    letters_a = analyze_string(chemical_a)
    letters_b = analyze_string(chemical_b)
    elements_a = get_number_chemical_elements(letters_a)
    elements_b = get_number_chemical_elements(letters_b)
    common_elements = get_number_common_elements(letters_a, letters_b)
    return round((common_elements / (elements_a + elements_b - common_elements)), 2)


def get_pivots(chemicals_length, number_processors):
    """Calculates the pivots to divide the chemicals between the threads."""
    pivots_list = []
    pivots_list.append(0)
    for i in range(number_processors - 1, 0, -1):
        pivots_list.append(
            int(round(chemicals_length - (math.sqrt(i / number_processors) * chemicals_length))))
    pivots_list.append(chemicals_length - 1)
    return pivots_list


def fill_compared_list(chemicals_list, pivot_min, pivot_max, compared_chemicals_list):
    """Fills a list with the comparison of the chemical compounds between two pivots in the list."""
    for i in range(pivot_min, pivot_max):
        for j in range(i + 1, len(chemicals_list)):
            coef = get_jac_tan_coefficient(chemicals_list[i][1], chemicals_list[j][1])
            row = chemicals_list[i][0] + "\t" + chemicals_list[j][0] + "\t" + str(coef)
            compared_chemicals_list.append(row)


def start_join_all(thread_list):
    """Starts and joins all the threads in the list."""
    for thread in thread_list:
        thread.start()
    for thread in thread_list:
        thread.join()


def write_file(number_threads, compared_chemicals_list, total_time):
    """Writes a tsv file with the compared chemicals."""
    with open("chem_sim_total_Python_Threads.tsv", "w") as record_file:
        record_file.write("Chem_ID_1\tChem_ID_2\tTanimoto_similarity\n")
        for i in range(number_threads):
            print(len(compared_chemicals_list[i]))
            for chemical in compared_chemicals_list[i]:
                record_file.write(chemical + "\n")
        record_file.write("Total time = " + str(total_time) + " [s]\n")


def print_to_console(number_threads, compared_chemicals_list, total_time):
    """Prints in console the compared chemicals."""
    counter = 0
    for i in range(number_threads):
        for chemical in compared_chemicals_list[i]:
            counter += 1
            print(chemical)
    print("Total Elements = " + str(counter))
    print("Total time = " + str(total_time) + " [s]")


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
    TOTAL_TIME = END_TIME - START_TIME
    write_file(NUMBER_THREADS, compared_chemicals, TOTAL_TIME)
    #print_to_console(NUMBER_THREADS, compared_chemicals, TOTAL_TIME)
