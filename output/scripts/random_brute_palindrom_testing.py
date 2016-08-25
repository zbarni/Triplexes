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
PATH_TEST_RANDOM = os.environ.get("TRIPLEXATOR_HOME") + "/output/test/results/random_5/"

sys.path.append(os.environ.get("TRIPLEXATOR_HOME") + '/python_bindings/')

import triplexator
import compare_brute_local_palindrom as compare
import time

# test variables
TESTS = 2
TEST_SIZE = 50000

# options
ERROR_RATE_LO = 10  # percentage
ERROR_RATE_HI = 20  # percentage
CONS_ERROR_LO = 1
CONS_ERROR_HI = 1
MIN_LENGTH_LO = 15
MIN_LENGTH_HI = 15
MAX_LENGTH    = ""

if __name__ == "__main__":
    with open(FILE_DNA_DATA_CHR1) as f:
        lines = f.readlines()
        genome = lines[1]

    for e in range(ERROR_RATE_LO, ERROR_RATE_HI + 1, 5):
        for c in range(CONS_ERROR_LO, CONS_ERROR_HI + 1):
            for l in range(MIN_LENGTH_LO, MIN_LENGTH_HI + 1, 5):

                # for local debug batch
                print("Doing batch with params: -e " + str(e) + " -c " + str(c) + " -l " + str(l) + MAX_LENGTH)
                res_file = open(PATH_TEST_RANDOM + 'random_l' + str(l) + '_e' + str(e) + '_c' + str(c) + '.log', 'w')
                success = True
                for test_nr in range(0, TESTS):
                    pos = randint(0, len(genome) - TEST_SIZE)
                    seq = genome[pos:pos + TEST_SIZE]

                    print("Running test #" + str(test_nr) + " ...", file=res_file)

                    # create new random data files
                    data_file_name = PATH_TEST_DATA + 'ss_' + str(test_nr) + '.data'
                    data_file = open(data_file_name, 'w')
                    print('>chr_test', file=data_file)
                    print(seq, file=data_file)
                    data_file.close()

                    print("Running bit-parallel local #" + str(test_nr) + " ...", file=res_file)
                    out_bit_file_name = "random_" + str(test_nr) + ".tpx"
                    dbg_bit_file_name = "random_" + str(test_nr) + ".dbg"
                    start_time = time.time()
                    triplexator.runTriplexator(
                        '-ss ' + data_file_name + " -ds " + data_file_name + " -fr off -fm 1 --bit-parallel-local " +
                        " -e " + str(e) + " -c " + str(c) + " -l " + str(l) + MAX_LENGTH +
                        " -od " + PATH_TEST_BPL + " -o " + out_bit_file_name + " &> " + dbg_bit_file_name)
                    print("BPL time: " + str(time.time() - start_time), file=res_file)

                    print("Running brute-force #" + str(test_nr) + " ...", file=res_file)
                    out_brute_file_name = "random_" + str(test_nr) + ".tpx"
                    dbg_brute_file_name = "random_" + str(test_nr) + ".dbg"
                    start_time = time.time()
                    triplexator.runTriplexator(
                        '-ss ' + data_file_name + " -ds " + data_file_name + " -fr off " +
                        " -e " + str(e) + " -c " + str(c) + " -l " + str(l) + MAX_LENGTH +
                        " -od " + PATH_TEST_BRUTE + " -o " + out_brute_file_name + " &> " + dbg_brute_file_name)
                    print("Brute time: " + str(time.time() - start_time), file=res_file)

                    (comp_res, comp_err) = compare.compare(PATH_TEST_BPL + out_bit_file_name, PATH_TEST_BRUTE + out_brute_file_name)
                    if not comp_res:
                        success = False
                        print("ERROR: " + comp_err, file=res_file)
                        exit(-1)
                if success:
                    print("SUCCESS", file=res_file)

