#!/usr/bin/python

"""
Main module to run various tests / benchmarks or submit jobs on the cluster, all related to Triplexator.
"""

from __future__ import print_function
from optparse import OptionParser
import os
import subprocess
import shlex
import data_analysis as da
import utils
import unit_tests as unit


def check_setup():
    """ Exit if TRIPLEXATOR_HOME environment variable is not set. """
    if os.environ.get("TRIPLEXATOR_HOME") is None:
        print("Please set TRIPLEXATOR_HOME environment variable properly!")
        print("Exiting..")
        exit(-1)


def get_triplexator_option(mode):
    """

    :param mode:
    :return:
    """
    if mode == "ab":
        return " "
    elif mode == "br":
        return ""
    elif mode == "bp":
        return "--bit-parallel"


def check_input_common(options, parser):
    # check mode
    if options.mode is None or (options.mode != "br" and options.mode != "bp" and options.mode != "ab"):
        parser.error("Mode is required / wrong mode given!")

    # if options is other then auto-binding, both TTS and TFO are required
    if options.mode == "br" or options.mode == 'bp':
        if options.inputTTS is None or options.inputTFO is None:
            parser.error("Both TTS and TFO input files are required!")
        if os.path.isfile(options.inputTTS) is False:
            parser.error("TTS file not found!")
        if os.path.isfile(options.inputTFO) is False:
            parser.error("TFO file not found!")

    # if max length is not given, set it to 30 (triplexator default)
    if options.maxLength is None:
        options.maxLength = 30

    # print some warning messages if some arguments are not given, and set them to default values
    if options.minLengthLow is None:
        print("WARNING: minimum triplex length not given, defaulting to " + utils.defaultMinLength)
        options.minLengthLow = utils.defaultMinLength
    if options.consLow is None:
        print("WARNING: lower bound on consecutive errors not given, defaulting to " + utils.defaultConsError)
        options.consLow = utils.defaultConsError
    if options.errLow is None:
        print("WARNING: lower bound on error rate not given, defaulting to " + utils.defaultErrorRate)
        options.errLow = utils.defaultErrorRate

    # if upper bound on some parameters is not given, set them to the corresponding lower bound
    if options.minLengthHigh is None:
        options.minLengthHigh = options.minLengthLow
    if options.consHigh is None:
        options.consHigh = options.consLow
    if options.errHigh is None:
        options.errHigh = options.errLow
    return options


def check_input_local(options, parser):
    """
    Sanitize and check input; make sure everything required is really there.
    :param options:
    :param parser: parser with read-in parameters
    :return: updated options
    """
    options.tpxOptions = " " + options.tpxOptions + " "
    options = check_input_common(options, parser)
    return options


def check_input_cluster(options, parser):
    """
    Sanitize and check input for the cluster; make sure everything required is really there.
    :param options:
    :param parser: parser with read-in parameters
    :return: updated options
    """
    options = check_input_common(options, parser)

    if options.clusterHour is None:
        parser.error("Job time limit is required!")
    if options.clusterMemory is None:
        parser.error("Job memory limit is required!")

    # check mode
    if options.mode is None or (options.mode != "br" and options.mode != "bp" and options.mode != "ab"):
        parser.error("Mode is required / wrong mode given!")

    return options


def run_cluster(options):
    """
    Iterates over the given parameter ranges and submits a job to the cluster for each parameter combination.
    :param options:
    :return:
    """
    # create path if non-existing
    if not os.path.isdir(options.dataOutDir):
        os.mkdir(options.dataOutDir)

    # set input files according to mode
    if options.mode == "br" or options.mode == "bp":
        input_file_options = '-ss ' + options.inputTFO + " -ds " + options.inputTTS + " "
    else:
        input_file_options = '-as ' + options.inputAB + " "

    for l in range(int(options.minLengthLow), int(options.minLengthHigh) + 1, 5):
        for e in range(int(options.errLow), int(options.errHigh) + 1, 5):
            for c in range(int(options.consLow), int(options.consHigh) + 1, 1):
                out_file_tpx = options.dataPrefix + '_l' + str(l) + '--' + str(options.maxLength) + '_e' + str(e)\
                                    + '_c' + str(c) + ".tpx"

                cmd = (utils.PATH_CLUSTER + "meta_bsub.sh bsub " + options.dataPrefix +
                       " " + options.dataOutDir + "/" + options.dataPrefix + '_l' + str(l) + '--' +
                       str(options.maxLength) + '_e' + str(e) + '_c' + str(c) + " " + options.clusterHour + " " +
                       options.clusterMemory + " " + utils.PATH_TRIPLEXATOR_HOME + " " + input_file_options + " " +
                       get_triplexator_option(options.mode) + " -e " + str(e) +
                       " -c " + str(c) + " -l " + str(l) + " -L " + options.maxLength + " " + " -od " +
                       options.dataOutDir + " -o " + out_file_tpx + options.tpxOptions)
                cmd = shlex.split(cmd)
                subprocess.call(cmd)


def run_locally(options):
    """
    Run triplexator locally by iterating over the given parameter ranges and submits a job to the cluster for each
    parameter combination.
    :param options:
    :return:
    """
    # create path if non-existing
    if not os.path.isdir(options.dataOutDir):
        os.mkdir(options.dataOutDir)

    # set input files according to mode
    if options.mode == "br" or options.mode == "bp":
        input_file_options = '-ss ' + options.inputTFO + " -ds " + options.inputTTS + " "
    else:
        input_file_options = '-as ' + options.inputAB + " "

    for l in range(int(options.minLengthLow), int(options.minLengthHigh) + 1, 5):
        for e in range(int(options.errLow), int(options.errHigh) + 1, 5):
            for c in range(int(options.consLow), int(options.consHigh) + 1, 1):
                out_file_tpx = utils.get_output_filename_from_parameters(options.dataPrefix, options.maxLength, l, e, c)

                utils.lib_triplexator.runTriplexator(input_file_options + get_triplexator_option(options.mode) +
                                                     " -e " + str(e) + " -c " + str(c) + " -l " + str(l) + " -L " +
                                                     str(options.maxLength) + " " + " -od " + options.dataOutDir +
                                                     " -o " + out_file_tpx + options.tpxOptions)


def create_parser():
    """
    Just adds all the options.
    :return: parser with options
    """
    parser = OptionParser()
    parser.add_option("-t", "--type", dest="type", help="where to run stuff", metavar=('local', 'cluster'))
    parser.add_option("--tts", dest="inputTTS", help="double stranded file (tts), usually dna", metavar="TTS")
    parser.add_option("--tfo", dest="inputTFO", help="single stranded file (tfo), usually rna", metavar="TFO")
    parser.add_option("--auto-binding-file", dest="inputAB",
                      help="single file (dna or rna) which will be searched against itself for triplexes", metavar="AB")
    parser.add_option("--error-low", dest="errLow", help="lower bound on error rate (%) [as a range]")
    parser.add_option("--error-high", dest="errHigh", help="upper bound on error rate (%) [as a range]")
    parser.add_option("--cons-low", dest="consLow", help="lower bound on #consecutive errors [as a range]")
    parser.add_option("--cons-high", dest="consHigh", help="upper bound on #consecutive errors [as a range]")
    parser.add_option("--min-length-low", dest="minLengthLow", help="lower bound on minimum length [as a range]")
    parser.add_option("--min-length-high", dest="minLengthHigh", help="upper bound on minimum length [as a range]")
    parser.add_option("--max-length", dest="maxLength", default="30", help="maximum triplex length")
    parser.add_option("--tpx-options", dest="tpxOptions", default="",
                      help="other options for triplexator, manually specified")
    parser.add_option("-m", "--mode", dest="mode", help="brute (br), bit-parallel (bp) or auto-binding (ab)",
                      metavar=('br', 'bp', 'ab'))

    parser.add_option("--data-prefix", dest="dataPrefix", default="",
                      help="prefix of input data files, e.g., files to be scanned in a directory")
    parser.add_option("--data-extension", dest="dataExtension", default="",
                      help="extension of input data files, e.g., files to be scanned in a directory")
    parser.add_option("--data-input-dir", dest="dataInDir", default="",
                      help="directory of input files, different from tfo and tts")
    parser.add_option("--data-output-dir", dest="dataOutDir", default="",
                      help="directory of output files, e.g., plots, etc.")

    parser.add_option("--data-analysis", dest="dataAnalysis", default=None,
                      metavar=('plot', 'len-dist', 'len-dist-multi'))
    parser.add_option("--length-distribution-triplexes", dest="lenDistTpx", default=None)
    parser.add_option("--length-distribution-regions", dest="lenDistRegions", default=None)
    parser.add_option("--jaccard-reference", dest="jaccardReference", default=None, metavar='jref')
    parser.add_option("--jaccard-query", dest="jaccardQuery", default=None, metavar='jquery')
    parser.add_option("--convert", dest="convert", default=None, metavar='tpx-to-bed')

    parser.add_option("--unit-tests", dest="unitTests", default=None, metavar=("test_directory_name", ""))

    parser.add_option("--cluster-time", dest="clusterHour", help="time for submitted job (cluster only)",
                      metavar="2:00")
    parser.add_option("--cluster-mem", dest="clusterMemory", help="memory for submitted job (cluster only)",
                      metavar="16000")
    parser.add_option("--valgrind", dest="valgrind", default=False, action="store_true", help="profile with valgrind",
                      metavar="TTS")
    return parser


if __name__ == "__main__":
    # import required stuff which wasn't possible at the beginning
    utils.lazy_imports()
    # some initial checking
    check_setup()

    m_parser = create_parser()
    (m_options, args) = m_parser.parse_args()

    if m_options.unitTests is not None:
        unit.run_bit_parallel_unit_test(m_options.unitTests)
    elif m_options.convert == "tpx-to-bed":
        utils.convert(m_options)
    elif m_options.dataAnalysis is not None:
        # unit.unit_data_analysis()
        da.main(m_options)
    else:
        if m_options.type == "local":
            m_options = check_input_local(m_options, m_parser)
            run_locally(m_options)
        elif m_options.type == "cluster":
            run_cluster(m_options)
