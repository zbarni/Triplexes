#!/usr/bin/python

"""

"""

from __future__ import print_function
from optparse import OptionParser
import sys
import os
import subprocess
import shlex
import data_analysis as da
import utils
import unit_tests as unit

result_dir = {"bpl": utils.PATH_DIR_BPL, "bp": utils.PATH_DIR_BP, "brute": utils.PATH_DIR_BRUTE}
# import compare_brute_local_palindrom as compare


def check_setup():
    if os.environ.get("TRIPLEXATOR_HOME") is None:
        print("Please set TRIPLEXATOR_HOME environment variable properly!")
        print("Exiting..")
        exit(-1)


def get_triplexator_option(mode):
    if mode == "bpl":
        return "--bit-parallel-local"
    elif mode == "brute":
        return ""
    elif mode == "bp":
        return "--bit-parallel"


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
    global result_dir

    out_dir = utils.PATH_CLUSTER + result_dir[options.mode]

    # create path if non-existing
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    for l in range(int(options.minLengthLow), int(options.minLengthHigh) + 1, 5):
        for e in range(int(options.errLow), int(options.errHigh) + 1, 5):
            for c in range(int(options.consLow), int(options.consHigh) + 1, 1):
                out_file_template = options.dataPrefix + '_l' + str(l) + '--' + str(options.maxLength) + '_e' + str(e) + '_c' + str(c)
                out_file_tpx = out_file_template + ".tpx"

                cmd = utils.PATH_CLUSTER + "meta_bsub.sh bsub " + options.dataPrefix + \
                      " " + out_dir + "/" + options.dataPrefix + '_l' + str(l) + '--' + str(options.maxLength) + '_e' +\
                      str(e) + '_c' + str(c) + " " + options.clusterHour + " " + options.clusterMemory + " " + \
                      utils.PATH_TRIPLEXATOR_HOME + " " + '-ss ' + options.inputTFO + " -ds " + options.inputTTS + " " +\
                      get_triplexator_option(options.mode) + " -e " + str(e) + " -c " + str(c) + " -l " + str(l) + \
                      " -L " + options.maxLength + " " + " -od " + utils.PATH_CLUSTER + result_dir[options.mode] + " -o " \
                      + out_file_tpx + options.tpxOptions
                cmd = shlex.split(cmd)
                subprocess.call(cmd)


def run_anything(options):
    # create path if non-existing
    if not os.path.isdir(options.dataOutDir):
        os.mkdir(options.dataOutDir)

    for l in range(int(options.minLengthLow), int(options.minLengthHigh) + 1, 5):
        for e in range(int(options.errLow), int(options.errHigh) + 1, 5):
            for c in range(int(options.consLow), int(options.consHigh) + 1, 1):
                out_file_tpx = get_output_filename(options, l, e, c)

                triplexator.runTriplexator('-ss ' + options.inputTFO + " -ds " + options.inputTTS + " " + \
                                           get_triplexator_option(options.mode) + " -e " + str(e) + " -c " + str(c) + \
                                           " -l " + str(l) + " -L " + options.maxLength + " " + " -od " + options.dataOutDir + \
                                           + " -o " + out_file_tpx + options.tpxOptions)


def run_test(options):
    global result_dir

    out_dir = utils.PATH_TEST + result_dir[options.mode]
    # create path if non-existing
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    for l in range(int(options.minLengthLow), int(options.minLengthHigh) + 1, 5):
        for e in range(int(options.errLow), int(options.errHigh) + 1, 5):
            for c in range(int(options.consLow), int(options.consHigh) + 1, 1):
                out_file_tpx = get_output_filename(options, l, e, c)

                if not options.valgrind:
                    triplexator.runTriplexator('-ss ' + options.inputTFO + " -ds " + options.inputTTS + " " + \
                                               get_triplexator_option(options.mode) + " -e " + str(e) + " -c " + str(c) + \
                                               " -l " + str(l) + " -L " + options.maxLength + " " + " -od " + utils.PATH_TEST + \
                                               result_dir[options.mode] + " -o " + out_file_tpx + options.tpxOptions)
                else:
                    cmd = "valgrind --tool=memcheck --leak-check=full " + utils.TRIPLEXATOR_BINARY + ' -ss ' +\
                          options.inputTFO + " -ds " + options.inputTTS + " " + \
                          get_triplexator_option(options.mode) + " -e " + str(e) + " -c " + str(c) + \
                          " -l " + str(l) + " -L " + options.maxLength + " " + " -od " + utils.PATH_TEST + \
                          result_dir[options.mode] + " -o " + out_file_tpx + options.tpxOptions
                    out_file_valgrind = get_valgrind_filename(options, l, e, c)
                    dump_file = open(out_dir + "/" + out_file_valgrind, 'w')
                    cmd = shlex.split(cmd)
                    subprocess.call(cmd, stdout=dump_file, stderr=dump_file)
                    dump_file.close()


def check_input(options, parser):
    """ """
    options.tpxOptions = " " + options.tpxOptions + " "
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
    """ Just adds all the options """
    parser = OptionParser()
    parser.add_option("-t", "--type", dest="type", help="what to run", metavar=('benchmark', 'cluster', 'test'))
    parser.add_option("--tts", dest="inputTTS", help="double stranded file (tts), usually dna", metavar="TTS")
    parser.add_option("--tfo", dest="inputTFO", help="single stranded file (tfo), usually rna", metavar="TFO")
    parser.add_option("--error-low", dest="errLow", help="lower bound on error rate (%) [as a range]")
    parser.add_option("--error-high", dest="errHigh", help="upper bound on error rate (%) [as a range]")
    parser.add_option("--cons-low", dest="consLow", help="lower bound on #consecutive errors [as a range]")
    parser.add_option("--cons-high", dest="consHigh", help="upper bound on #consecutive errors [as a range]")
    parser.add_option("--min-length-low", dest="minLengthLow", help="lower bound on minimum length [as a range]")
    parser.add_option("--min-length-high", dest="minLengthHigh", help="upper bound on minimum length [as a range]")
    parser.add_option("--max-length", dest="maxLength", default="30", help="maximum triplex length")
    parser.add_option("--tpx-options", dest="tpxOptions", default="", help="other options for triplexator, manually specified")
    parser.add_option("-m", "--mode", dest="mode", help="brute (), bit-parallel (bp) or bit-parallel-local (bpl)", metavar=('brute', 'bp', 'bpl'))

    parser.add_option("--data-prefix", dest="dataPrefix", default="", help="prefix of input data files, e.g., files to be scanned in a directory")
    parser.add_option("--data-extension", dest="dataExtension", default="", help="extension of input data files, e.g., files to be scanned in a directory")
    parser.add_option("--data-input-dir", dest="dataInDir", default="", help="directory of input files, different from tfo and tts")
    parser.add_option("--data-output-dir", dest="dataOutDir", default="", help="directory of output files, e.g., plots, etc.")

    parser.add_option("--data-analysis", dest="dataAnalysis", default=None, metavar=('plot'))
    parser.add_option("--convert", dest="convert", default=None, metavar=('tpx-to-bed'))

    parser.add_option("--unit-tests", dest="unitTests", default=None, metavar=("test_directory_name", ""))

    parser.add_option("--cluster-time", dest="clusterHour", help="time for submitted job (cluster only)", metavar="2:00")
    parser.add_option("--cluster-mem", dest="clusterMemory", help="memory for submitted job (cluster only)", metavar="16000")
    parser.add_option("--valgrind", dest="valgrind", default=False, action="store_true", help="profile with valgrind", metavar="TTS")
    return parser

if __name__ == "__main__":
    utils.lazy_imports()
    check_setup()
    m_parser = create_parser()
    (m_options, args) = m_parser.parse_args()

    if m_options.unitTests is not None:
        unit.run_bit_parallel_unit_test(m_options.unitTests)
    elif m_options.convert == "tpx-to-bed":
        utils.convert(m_options)
    elif m_options.dataAnalysis is not None:
        da.do_data_analysis(m_options)
    else:
        m_options = check_input(m_options, m_parser)
        if m_options.type == "benchmark":
            run_anything(m_options)
        elif m_options.type == "cluster":
            run_cluster(m_options)
        elif m_options.type == "test":
            run_test(m_options)
