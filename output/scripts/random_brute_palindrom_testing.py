#!/usr/bin/python

from __future__ import print_function
import sys, os
from random import randint

FILE_DNA_DATA_CHR1 = os.environ.get("TRIPLEXATOR_HOME") + "/data/dna/mm9/mm9.chr1.oneline.fa"
PATH_TEST_DATA = os.environ.get("TRIPLEXATOR_HOME") + "/output/test/data/"
PATH_TEST_OUTPUT = os.environ.get("TRIPLEXATOR_HOME") + "/output/test/output/"

sys.path.append(os.environ.get("TRIPLEXATOR_HOME") + 'python_bindings/')
import triplexator


TESTS = 2
TEST_SIZE = 50000

if __name__ == "__main__":
    with open(FILE_DNA_DATA_CHR1) as f:
        lines = f.readlines()
        genome = lines[1]

    for test_nr in range(0, TESTS):
        pos = randint(0, len(genome) - TEST_SIZE)
        seq = genome[pos:pos + TEST_SIZE]

        data_file = open(PATH_TEST_DATA + 'ss_' + str(test_nr) + '.data', 'w')
        print('>chr_test', file=data_file)
        print(seq, file=data_file)
        data_file.close()
        triplexator.runTriplexator('')
