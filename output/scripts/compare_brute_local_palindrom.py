#!/usr/bin/python

import sys
MAX_SHIFT_OFFSET = 0


def get_next_brute_line(idx, content):
    global MAX_SHIFT_OFFSET
    while True and idx < len(content):
        line = content[idx].split('\t')
        if line[0][0] == '#':
            return None

        tfo_beg = int(line[1])
        tfo_end = int(line[2])

        tts_beg = int(line[4])
        tts_end = int(line[5])

        length = tts_end - tts_beg

        if abs(tts_end - tfo_end) <= length + MAX_SHIFT_OFFSET:
            return idx
        idx += 1

    return None


def matching_palindrom(idx, content, reference_line):
    if idx >= len(content):
        return False
    print "comparing..."
    print reference_line
    print content[idx]
    print "----"
    return reference_line == content[idx]

def compare(pal_filename, bru_filename):
    with open(pal_filename) as f:
        pal_content = f.readlines()

    with open(bru_filename) as f:
        bru_content = f.readlines()

    # init indices and skip first line (not needed)
    pal_idx = bru_idx = 0
    while pal_idx < len(pal_content) and bru_idx < len(bru_content):
        pal_line = pal_content[pal_idx]
        bru_idx = get_next_brute_line(bru_idx, bru_content)
        if bru_idx is None:
            break
        # print "found one: ", bru_content[bru_idx]
        if not matching_palindrom(pal_idx, pal_content, bru_content[bru_idx]):
            print "ERROR - not matching palindrom found for line below:"
            print "\tbrute line: ", bru_content[bru_idx]
            exit(-1)

        bru_idx += 1
        pal_idx += 1

    print ">>> SUCCESS."


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "This script requires 2 txp files as parameters, local_myers and brute_force outputs."
        exit(-1)

    compare(sys.argv[1], sys.argv[2])
