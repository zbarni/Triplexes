#!/usr/bin/python
import numpy as np
import os
import re
from matplotlib import pyplot as plt
import pickle
import pprint
# import utils

pp = pprint.PrettyPrinter(indent=2)


def plot_length_distribution(tpx_content, filename, options):
    lengths = {}

    for line in tpx_content:
        cols = line.split('\t')
        if line[0] == '#':
            continue

        tfo_start = int(cols[1])
        tfo_end = int(cols[2])
        tfo_len = tfo_end - tfo_start
        lengths[tfo_len] = 1 if tfo_len not in lengths else lengths[tfo_len] + 1

    out_filename = filename + ".len.serialized"

    if not os.path.isdir(options.dataOutDir):
        os.mkdir(options.dataOutDir)

    pickle.dump(lengths, open(options.dataOutDir + "/" + out_filename, "wb"))


def do_data_analysis(options):
    for f in os.listdir(options.dataInDir):
        if f.startswith(options.dataPrefix) and f.endswith(".tpx"):
            full_file_path = options.dataInDir + '/' + f
            if os.path.isfile(full_file_path):
                tpx_content = (open(full_file_path)).readlines()
                plot_length_distribution(tpx_content, f, options)
