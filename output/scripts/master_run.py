#!/usr/bin/python

"""

"""

from __future__ import print_function
import sys, os
from random import randint
from optparse import OptionParser
import subprocess
import shlex
import data_analysis as da

if os.environ.get("TRIPLEXATOR_HOME") is None:
    print("Please set TRIPLEXATOR_HOME environment variable properly!")
    print("Exiting..")
    exit(-1)

FILE_DNA_DATA_CHR1 = os.environ.get("TRIPLEXATOR_HOME") + "/data/dna/mm9/mm9.chr1.oneline.fa"
PATH_TRIPLEXATOR_HOME = os.environ.get("TRIPLEXATOR_HOME")
TRIPLEXATOR_BINARY = PATH_TRIPLEXATOR_HOME + "/triplexator/bin/triplexator"
PATH_TEST_DATA = os.environ.get("TRIPLEXATOR_HOME") + "/output/test/data/"
PATH_TEST_BPL = os.environ.get("TRIPLEXATOR_HOME") + "/output/test/results/bit_parallel_local/"
PATH_TEST_BP = os.environ.get("TRIPLEXATOR_HOME") + "/output/test/results/"
PATH_TEST_BRUTE = os.environ.get("TRIPLEXATOR_HOME") + "/output/test/results/brute/"
PATH_TEST_RANDOM = os.environ.get("TRIPLEXATOR_HOME") + "/output/test/results/random_4/"
PATH_CLUSTER = os.environ.get("TRIPLEXATOR_HOME") + "/output/cluster/"
PATH_TEST = os.environ.get("TRIPLEXATOR_HOME") + "/output/test/"
PATH_DIR_BPL = "bpl/"
PATH_DIR_BRUTE = "brute/"
PATH_DIR_BP = "bp/"

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

result_dir = {"bpl": PATH_DIR_BPL, "bp": PATH_DIR_BP, "brute": PATH_DIR_BRUTE}


def get_triplexator_option(mode):
    if mode == "bpl":
        return "--bit-parallel-local"


def get_output_filename(options, l, e, c):
    out_file_template = "valgrind_" if options.valgrind else ""
    out_file_template += options.dataPrefix + '_l' + str(l) + '--' + str(options.maxLength) + '_e' + str(e) + '_c' + str(c)
    out_file_tpx = out_file_template + ".tpx"
    return out_file_tpx


def get_valgrind_filename(options, l, e, c):
    out_file_template = "valgrind_" if options.valgrind else ""
    out_file_template += options.dataPrefix + '_l' + str(l) + '--' + str(options.maxLength) + '_e' + str(e) + '_c' + str(c)
    out_file_val = out_file_template + ".val"
    return out_file_val


def run_cluster(options):
    global PATH_CLUSTER
    global PATH_TRIPLEXATOR_HOME
    global result_dir

    out_dir = PATH_CLUSTER + result_dir[options.mode]

    # create path if non-existing
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    for l in range(int(options.minLengthLow), int(options.minLengthHigh) + 1, 5):
        for e in range(int(options.errLow), int(options.errHigh) + 1, 5):
            for c in range(int(options.consLow), int(options.consHigh) + 1, 1):
                out_file_template = options.dataPrefix + '_l' + str(l) + '--' + str(options.maxLength) + '_e' + str(e) + '_c' + str(c)
                out_file_log = out_file_template + ".log"
                out_file_tpx = out_file_template + ".tpx"

                cmd = PATH_CLUSTER + "meta_bsub.sh bsub " + options.dataPrefix + \
                      " " + out_dir + "/" + options.dataPrefix + '_l' + str(l) + '--' + str(options.maxLength) + '_e' + str(e) + '_c' + str(c) + \
                      " " + options.clusterHour + " " + options.clusterMemory + " " +\
                      PATH_TRIPLEXATOR_HOME + " " + '-ss ' + options.inputTFO + " -ds " + options.inputTTS + " " + \
                      get_triplexator_option(options.mode) + " -e " + str(e) + " -c " + str(c) + " -l " + str(l) + \
                      " -L " + options.maxLength + " " + " -od " + PATH_CLUSTER + result_dir[options.mode] + " -o " \
                      + out_file_tpx
                cmd = shlex.split(cmd)
                subprocess.call(cmd)


def run_test(options):
    global PATH_TEST
    global PATH_TRIPLEXATOR_HOME
    global TRIPLEXATOR_BINARY
    global result_dir

    out_dir = PATH_TEST + result_dir[options.mode]

    # create path if non-existing
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    for l in range(int(options.minLengthLow), int(options.minLengthHigh) + 1, 5):
        for e in range(int(options.errLow), int(options.errHigh) + 1, 5):
            for c in range(int(options.consLow), int(options.consHigh) + 1, 1):
                # out_file_log = out_file_template + ".log"
                out_file_tpx = get_output_filename(options, l, e, c)

                if not options.valgrind:
                    triplexator.runTriplexator('-ss ' + options.inputTFO + " -ds " + options.inputTTS + " " + \
                                               get_triplexator_option(options.mode) + " -e " + str(e) + " -c " + str(c) + \
                                               " -l " + str(l) + " -L " + options.maxLength + " " + " -od " + PATH_TEST + \
                                               result_dir[options.mode] + " -o " + out_file_tpx)
                else:
                    cmd = "valgrind --tool=memcheck --leak-check=full " + TRIPLEXATOR_BINARY + ' -ss ' +\
                          options.inputTFO + " -ds " + options.inputTTS + " " + \
                          get_triplexator_option(options.mode) + " -e " + str(e) + " -c " + str(c) + \
                          " -l " + str(l) + " -L " + options.maxLength + " " + " -od " + PATH_TEST + \
                          result_dir[options.mode] + " -o " + out_file_tpx
                    out_file_valgrind = get_valgrind_filename(options, l, e, c)
                    dump_file = open(out_dir + "/" + out_file_valgrind, 'w')
                    cmd = shlex.split(cmd)
                    subprocess.call(cmd, stdout=dump_file, stderr=dump_file)
                    dump_file.close()


def check_input(options, parser):
    # checks for data analysis
    if options.dataAnalysis:
        return options

    ####################################################################
    if options.inputTTS is None or options.inputTFO is None:
        parser.error("Both TTS and TFO input files are required!")
    if os.path.isfile(options.inputTTS) is False:
        parser.error("TTS file not found!")
    if os.path.isfile(options.inputTFO) is False:
        parser.error("TFO file not found!")

    if options.mode is None or (options.mode != "brute" and options.mode != "bp" and options.mode != "bpl"):
        parser.error("Mode is required / wrong mode given!")
        return

    if options.minLengthHigh is None:
        options.minLengthHigh = options.minLengthLow

    if options.maxLength is None:
        options.maxLength = 30
    return options

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
    parser.add_option("--max-length", dest="maxLength", default="30", help="double stranded file (tts), usually dna", metavar="TTS")
    parser.add_option("-m", "--mode", dest="mode", help="double stranded file (tts), usually dna", metavar="TTS")
    parser.add_option("--data-prefix", dest="dataPrefix", default="", help="double stranded file (tts), usually dna", metavar="TTS")
    parser.add_option("--data-input-dir", dest="dataInDir", default="", help="double stranded file (tts), usually dna",
                      metavar="TTS")
    parser.add_option("--data-output-dir", dest="dataOutDir", default="", help="double stranded file (tts), usually dna")
    parser.add_option("--data-analysis", dest="dataAnalysis", default=False, action="store_true")
    parser.add_option("--cluster-time", dest="clusterHour", help="double stranded file (tts), usually dna", metavar="TTS")
    parser.add_option("--cluster-mem", dest="clusterMemory", help="double stranded file (tts), usually dna", metavar="TTS")
    parser.add_option("--valgrind", dest="valgrind", default=False, action="store_true", help="double stranded file (tts), usually dna", metavar="TTS")
    return parser

if __name__ == "__main__":
    parser = create_parser()
    (options, args) = parser.parse_args()

    options = check_input(options, parser)

    if options.dataAnalysis:
        da.do_data_analysis(options)

    elif options.type == "benchmark":
        pass
    elif options.type == "cluster":
        run_cluster(options)
    elif options.type == "test":
        run_test(options)
