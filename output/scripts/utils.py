import os
import gc
import sys

if os.environ.get("TRIPLEXATOR_HOME") is None:
    print("Please set TRIPLEXATOR_HOME environment variable properly!")
    print("Exiting..")
    exit(-1)

FILE_DNA_DATA_CHR1 = os.environ.get("TRIPLEXATOR_HOME") + "/data/dna/mm9/mm9.chr1.oneline.fa"
PATH_TRIPLEXATOR_HOME = os.environ.get("TRIPLEXATOR_HOME")
PATH_RGT_HOME = os.environ.get("RGT_HOME")
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

import_RNADNABindingSet = None

if PATH_RGT_HOME is not None:
    sys.path.append(PATH_RGT_HOME)
    import GenomicRegionSet as grs
    from triplex import RNADNABindingSet as _RNADNABindingSet
    import_RNADNABindingSet = _RNADNABindingSet


def get_files_with_extension(directory, prefix="", ext=""):
    if not os.path.isdir(directory):
        raise

    files = []
    for f in os.listdir(directory):
        if f.startswith(prefix) and f.endswith(ext):
            files.append(directory + '/' + f)
    return files


def convert_file(input_file_path, output_file_path):
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
    """ """
    if options.convert == "tpx-to-bed":
        files = get_files_with_extension(options.dataInDir, options.dataPrefix, options.dataExtension)
        for f in files:
            print("Converting file " + f)
            dbs_out_filename = f.replace(options.dataExtension, '.tts.bed')
            convert_file(f, dbs_out_filename)
            print("Goodbye Conversion.")
            # gc.collect()
    else:
        print("Invalid convert option.")

