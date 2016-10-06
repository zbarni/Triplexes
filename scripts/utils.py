"""
Utils module containing some global path variables and useful functions.
"""

import os
import sys

PATH_TRIPLEXATOR_HOME   = os.environ.get("TRIPLEXATOR_HOME")
PATH_RGT_HOME           = os.environ.get("RGT_HOME")
TRIPLEXATOR_BINARY      = PATH_TRIPLEXATOR_HOME + "/triplexator/bin/triplexator"
PATH_TEST_DATA          = os.environ.get("TRIPLEXATOR_HOME") + "/results/test/data/"
PATH_TEST_BPL           = os.environ.get("TRIPLEXATOR_HOME") + "/results/test/results/bit_parallel_local/"
PATH_TEST_BP            = os.environ.get("TRIPLEXATOR_HOME") + "/results/test/results/"
PATH_TEST_BRUTE         = os.environ.get("TRIPLEXATOR_HOME") + "/results/test/results/brute/"
PATH_TEST_RANDOM        = os.environ.get("TRIPLEXATOR_HOME") + "/results/test/results/random_4/"
PATH_CLUSTER            = os.environ.get("TRIPLEXATOR_HOME") + "/results/cluster/"
PATH_TEST               = os.environ.get("TRIPLEXATOR_HOME") + "/results/test/"
PATH_UNIT_TESTS         = os.environ.get("TRIPLEXATOR_HOME") + "/results/unit_tests/"
PATH_PYTHON_BINDINGS    = PATH_TRIPLEXATOR_HOME + '/python_bindings/'

defaultMinLength = 15
defaultConsError = 1
defaultErrorRate = 20

lib_RNADNABindingSet = None
lib_triplexator = None


def lazy_imports():
    """
    Do lazy import of triplexator library and some RGT library.
    :return:
    """
    global PATH_RGT_HOME
    global PATH_TRIPLEXATOR_HOME
    global PATH_PYTHON_BINDINGS

    # triplexator python binding
    if not os.path.isdir(PATH_PYTHON_BINDINGS):
        raise ImportError('Python bindings directory not found, make sure you checked out the repo correctly!')
    sys.path.append(PATH_TRIPLEXATOR_HOME + '/python_bindings/')
    import triplexator as _triplexator
    global lib_triplexator
    lib_triplexator = _triplexator


def get_output_filename_from_parameters(prefix, lmin, lmax, e, c):
    # out_file_template = "valgrind_" if options.valgrind else ""
    out_file_template = prefix + '_l' + str(lmin) + '--' + str(lmax) + '_c' + str(c) + '_e' + str(e)
    out_file_tpx = out_file_template + ".tpx"
    return out_file_tpx


def get_valgrind_filename(options, l, e, c):
    out_file_template = "valgrind_" if options.valgrind else ""
    out_file_template += options.dataPrefix + '_l' + str(l) + '--' + str(options.maxLength) + '_e' + str(e) + '_c' + str(c)
    out_file_val = out_file_template + ".val"
    return out_file_val


def get_files_with_extension(directory, prefix="", ext=""):
    if not os.path.isdir(directory):
        raise

    files = []
    for f in os.listdir(directory):
        if f.startswith(prefix) and f.endswith(ext):
            files.append(directory + '/' + f)
    return files


def get_tpx_lines(full_file_path):
    """

    :param full_file_path:
    :return:
    """
    result = []
    with open(full_file_path) as input_file:
        for line in input_file.readlines():
            if line[0] == "#":
                continue
            result.append(line)
    return result


def get_tpx_lines_by_chromosome(full_file_path):
    """

    :param full_file_path:
    :return:
    """
    result = {}
    with open(full_file_path) as input_file:
        for line in input_file.readlines():
            if line[0] == "#":
                continue
            cols = line.split('\t')
            if cols[0] not in result:
                result[cols[0]] = []
            result[cols[0]].append(line)
    return result


def clean_nonmatching_chromosomes(full_file_path):
    """
    Remove triplexes which are located on different chromosomes.
    :param full_file_path:
    :param tpx_content:
    :return:
    """
    print ("Cleaning " + full_file_path)

    tpx_content = get_tpx_lines(full_file_path)
    tpx_file_write = open(full_file_path, 'wb')

    cleaned_content = []
    for line in tpx_content:
        cols = line.split('\t')
        if line[0] != '#' and cols[0] != cols[3]:
            continue
        annotated_line = line.replace(cols[0], cols[0]+":0") if ':' not in cols[0] else line
        # annotated_line = line.replace(":0:0", ":0")  # if ':' not in cols[0] else line
        cleaned_content.append(annotated_line)

    tpx_file_write.writelines(cleaned_content)
    tpx_file_write.close()


def convert_file(input_file_path, output_file_path):
    """

    :param input_file_path:
    :param output_file_path:
    :return:
    """
    bed_file = open(output_file_path, 'w')
    with open(input_file_path) as input_file:
        for line in input_file.readlines():
            if line[0] == "#":
                continue
            cols = line.split('\t')
            (tfo_chr, tfo_start_offset) = cols[0].split(':')
            (tts_chr, tts_start_offset) = cols[3].split(':')
            (tfo_start_pos, tfo_end_pos, tts_start_pos, tts_end_pos) = (cols[1], cols[2], cols[4], cols[5])
            assert (tfo_chr == tts_chr)

            separator = "\t"
            parallel = cols[11]
            match_length = cols[6]
            strand = cols[10]

            out_line = tfo_chr + separator + str(int(tts_start_pos) + int(tts_start_offset)) + separator + str(int(tts_end_pos) + int(tts_start_offset))
            out_line += separator + tfo_start_pos + "-" + str(int(tfo_end_pos) + int(tfo_start_offset)) + "-" + parallel
            out_line += separator + match_length + separator + strand + "\n"

            bed_file.write(out_line)
    bed_file.close()


def convert(options):
    """

    :param options:
    :return:
    """
    if options.convert == "tpx-to-bed":
        files = get_files_with_extension(options.dataInDir, options.dataPrefix, options.dataExtension)

        for f in files:
            print("Converting file " + f)
            ext = options.dataExtension
            if options.dataExtension is None or options.dataExtension == "":
                ext = f.split('.')[-1]

            dbs_out_filename = f.replace(ext, 'tts.bed')
            # clean_nonmatching_chromosomes(f)
            convert_file(f, dbs_out_filename)
            print("Goodbye Conversion.")
    else:
        print("Invalid convert option.")


def check_io_files(files):
    if files is None or type(files) is not list:
        raise IOError("Invalid files parameter!")

    for f in files:
        if f is None:
            raise IOError("Invalid files parameter!")
        if ' ' in f:
            for ff in f.split(' '):
                if not os.path.isfile(ff):
                    raise IOError("File " + ff + " doesn't exist!")
        else:
            if not os.path.isfile(f):
                raise IOError("File " + f + " doesn't exist!")