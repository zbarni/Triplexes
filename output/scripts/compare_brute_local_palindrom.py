#!/usr/bin/python

"""
This module compares the outputs of 2 .tpx files (usually bit-parallel-local and brute-force) and checks if the
  all the semi-palindromic regions found by the original brute-force version are contained in the bit-parallel output.

Author: Barna Zajzon
"""

import sys, os
from subprocess import call
MAX_SHIFT_OFFSET = 0


def get_next_brute_line(idx, content):
    global MAX_SHIFT_OFFSET
    while True and idx < len(content):
        line = content[idx].split('\t')
        if line[0][0] == '#':
            return None

        tfo_end = int(line[2])

        tts_beg = int(line[4])
        tts_end = int(line[5])

        length = tts_end - tts_beg

        # print "get next: " + str(idx + 1) + ", val: " + str(abs(tts_end - tfo_end))
        if abs(tts_end - tfo_end) <= length + MAX_SHIFT_OFFSET:
            return idx
        idx += 1

    return None


def matching_palindrom(pal_idx, pal_content, reference_line):
    if pal_idx >= len(pal_content):
        return False
    return reference_line == pal_content[pal_idx]


def compare(pal_filename, bru_filename):
    """
    Compares the outputs of the two versions. If all matches of brute-force are contained in the bit-parallel-local as
     well, SUCCESS is printed. Otherwise program is exited and the line containing the error is printed.

    :param pal_filename: .tpx output of bit-parallel-local version
    :param bru_filename: .tpx output of brute-force version
    :return: nothing
    """
    tmp_pal_filename = "/tmp/bit-parallel-local"
    tmp_bru_filename = "/tmp/brute-force"

    # sort both versions
    with open(tmp_pal_filename, 'w') as out:
        call(["sort", pal_filename], stdout=out)
    with open(tmp_bru_filename, 'w') as out:
        call(["sort", bru_filename], stdout=out)

    # open and read files
    with open(tmp_pal_filename) as f:
        pal_content = f.readlines()

    with open(tmp_bru_filename) as f:
        bru_content = f.readlines()

    # delete temp files
    # os.remove(tmp_pal_filename)
    # os.remove(tmp_bru_filename)

    # init indices and skip first line (not needed)
    pal_idx = bru_idx = 0
    while pal_idx < len(pal_content) and bru_idx < len(bru_content):
        # get index of brute-force line where next match lies (then search for it in the palindromic)
        bru_idx = get_next_brute_line(bru_idx, bru_content)
        if bru_idx is None:
            break

        # print "> " + str(bru_idx + 1) + " vs " + str(pal_idx + 1)
        if not matching_palindrom(pal_idx, pal_content, bru_content[bru_idx]):
            print "ERROR - mismatching palindrom found for line below:"
            print "#" + str(bru_idx + 1) + "\tbrute line: ", bru_content[bru_idx]
            print "#" + str(pal_idx + 1) + "\tbit-parallel line: ", pal_content[pal_idx]
            print "!!! ERROR."
            return
            # exit(-1)

        bru_idx += 1
        pal_idx += 1

    print ">>> SUCCESS."


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "This script requires 2 txp files as parameters, local_myers and brute_force outputs."
        exit(-1)

    compare(sys.argv[1], sys.argv[2])
