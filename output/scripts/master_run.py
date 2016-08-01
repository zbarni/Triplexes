#!/usr/bin/python

"""

"""

from __future__ import print_function
import sys, os
from random import randint
from optparse import OptionParser
import subprocess

if os.environ.get("TRIPLEXATOR_HOME") is None:
    print("Please set TRIPLEXATOR_HOME environment variable properly!")
    print("Exiting..")
    exit(-1)

FILE_DNA_DATA_CHR1 = os.environ.get("TRIPLEXATOR_HOME") + "/data/dna/mm9/mm9.chr1.oneline.fa"
PATH_TEST_DATA = os.environ.get("TRIPLEXATOR_HOME") + "/output/test/data/"
PATH_TEST_BPL = os.environ.get("TRIPLEXATOR_HOME") + "/output/test/results/bit_parallel_local/"
PATH_TEST_BP = os.environ.get("TRIPLEXATOR_HOME") + "/output/test/results/"
PATH_TEST_BRUTE = os.environ.get("TRIPLEXATOR_HOME") + "/output/test/results/brute/"
PATH_TEST_RANDOM = os.environ.get("TRIPLEXATOR_HOME") + "/output/test/results/random_4/"
PATH_CLUSTER = os.environ.get("TRIPLEXATOR_HOME") + "/output/cluster/"
PATH_DIR_BPL = "palindrom/"
PATH_DIR_BRUTE = "brute/"
PATH_DIR_BP = "/"


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
CONS_ERROR_HI = 3
MIN_LENGTH_LO = 15
MIN_LENGTH_HI = 25
MAX_LENGTH    = ""


def run_cluster(options, parser):
    global PATH_CLUSTER
    global PATH_DIR_BP
    global PATH_DIR_BPL
    global PATH_DIR_BRUTE

    if options.mode is None or (options.mode != "brute" and options.mode != "bp" and options.mode != "bpl"):
        parser.error("Mode is required / wrong mode given!")
        return
    outDir = PATH_CLUSTER + options.mode


def check_input(options, parser):
    if options.inputTTS is None or options.inputTFO is None:
        parser.error("Both TTS and TFO input files are required!")
    if os.path.isfile(options.inputTTS) is False:
        parser.error("TTS file not found!")
    if os.path.isfile(options.inputTFO) is False:
        parser.error("TFO file not found!")


def create_parser():
    parser = OptionParser()
    parser.add_option("-t", "--type", dest="type", help="what to run", metavar=('benchmark', 'cluster'))
    parser.add_option("--tts", dest="inputTTS", help="double stranded file (tts), usually dna", metavar="TTS")
    parser.add_option("--tfo", dest="inputTFO", help="single stranded file (tfo), usually rna", metavar="TFO")
    parser.add_option("--error-low", dest="errLow", help="double stranded file (tts), usually dna", metavar="TTS")
    parser.add_option("--error-high", dest="errHigh", help="double stranded file (tts), usually dna", metavar="TTS")
    parser.add_option("--cons-low", dest="consLow", help="double stranded file (tts), usually dna", metavar="TTS")
    parser.add_option("--cons-high", dest="consHigh", help="double stranded file (tts), usually dna", metavar="TTS")
    parser.add_option("--min-length-low", dest="minLengthLow", help="double stranded file (tts), usually dna",
                      metavar="TTS")
    parser.add_option("--min-length-high", dest="minLengthHigh", help="double stranded file (tts), usually dna",
                      metavar="TTS")
    parser.add_option("--max-length", dest="maxLength", help="double stranded file (tts), usually dna", metavar="TTS")
    parser.add_option("-m", "--mode", dest="mode", help="double stranded file (tts), usually dna", metavar="TTS")
    return parser

if __name__ == "__main__":
    parser = create_parser()
    (options, args) = parser.parse_args()

    check_input(options, parser)


    if options.type == "benchmark":
        pass
    elif options.type == "cluster":
        run_cluster(options, parser)
