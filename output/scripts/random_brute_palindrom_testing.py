#!/usr/bin/python

"""

"""

from __future__ import print_function
import sys, os
from random import randint

if os.environ.get("TRIPLEXATOR_HOME") is None:
    print("Please set TRIPLEXATOR_HOME environment variable properly!")
    print("Exiting..")
    exit(-1)

FILE_DNA_DATA_CHR1 = os.environ.get("TRIPLEXATOR_HOME") + "/data/dna/mm9/mm9.chr1.oneline.fa"
PATH_TEST_DATA = os.environ.get("TRIPLEXATOR_HOME") + "/output/test/data/"
PATH_TEST_BPL = os.environ.get("TRIPLEXATOR_HOME") + "/output/test/results/bit_parallel_local/"
PATH_TEST_BP = os.environ.get("TRIPLEXATOR_HOME") + "/output/test/results/"
PATH_TEST_BRUTE = os.environ.get("TRIPLEXATOR_HOME") + "/output/test/results/brute/"

sys.path.append(os.environ.get("TRIPLEXATOR_HOME") + '/python_bindings/')

import triplexator
import compare_brute_local_palindrom as compare

# test variables
TESTS = 1
TEST_SIZE = 10000

# options
ERROR_RATE = " -e 20 "  # percentage
CONS_ERROR = " -c 2 "
MIN_LENGTH = " -l 15 "
MAX_LENGTH = " -L -1 "

if __name__ == "__main__":
    with open(FILE_DNA_DATA_CHR1) as f:
        lines = f.readlines()
        genome = lines[1]

    for test_nr in range(0, TESTS):
        pos = randint(0, len(genome) - TEST_SIZE)
        seq = genome[pos:pos + TEST_SIZE]

        print("Running test #" + str(test_nr) + " ...")

        # create new random data files
        data_file_name = PATH_TEST_DATA + 'ss_' + str(test_nr) + '.data'
        # data_file = open(data_file_name, 'w')
        # print('>chr_test', file=data_file)
        # print(seq, file=data_file)
        # data_file.close()

        print("Running bit-parallel local #" + str(test_nr) + " ...")
        out_bit_file_name = "random_" + str(test_nr) + ".tpx"
        dbg_bit_file_name = "random_" + str(test_nr) + ".dbg"
        triplexator.runTriplexator(
            '-ss ' + data_file_name + " -ds " + data_file_name + " -fr off -fm 1 --bit-parallel-local " +
            ERROR_RATE + CONS_ERROR + MIN_LENGTH + MAX_LENGTH +
            " -od " + PATH_TEST_BPL + " -o " + out_bit_file_name + " &> " + dbg_bit_file_name)

        print("Running brute-force #" + str(test_nr) + " ...")
        out_brute_file_name = "random_" + str(test_nr) + ".tpx"
        dbg_brute_file_name = "random_" + str(test_nr) + ".dbg"
        triplexator.runTriplexator(
            '-ss ' + data_file_name + " -ds " + data_file_name + " -fr off " +
            ERROR_RATE + CONS_ERROR + MIN_LENGTH + MAX_LENGTH +
            " -od " + PATH_TEST_BRUTE + " -o " + out_brute_file_name + " &> " + dbg_brute_file_name)
        compare.compare(PATH_TEST_BPL + out_bit_file_name, PATH_TEST_BRUTE + out_brute_file_name)
        print("-----\n")
