import os
import sys
import intervaltree
import random
# from matplotlib import pyplot as plt
import pprint
import utils
import time

lib_RNADNABindingSet = None
lib_viz = None

pp = pprint.PrettyPrinter(indent=2)


def lazy_import():
    # RGT library
    if utils.PATH_RGT_HOME is not None:
        sys.path.append(utils.PATH_RGT_HOME)
        global lib_RNADNABindingSet
        global lib_viz

        from triplex import RNADNABindingSet as _RNADNABindingSet
        import viz as lib_viz

        lib_RNADNABindingSet = _RNADNABindingSet


def _get_tpx_lines(full_file_path):
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


def _get_tpx_lines_by_chromosome(full_file_path):
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


# def plot_length_distributions(full_file_paths, filenames, options):
#     """
#
#     :param full_file_paths:
#     :param filenames:
#     :param options:
#     :return:
#     """
#     plt.clf()
#     # out_filename = filenames[0] + ".combined.plot.png"
#     out_filename = filenames[0] + ".plot.png"
#     markers = ["o", "s"]
#     colors = ["green", "red"]
#     for key, filename in enumerate(filenames):
#         full_file_path = full_file_paths[key]
#         print ("Plotting " + filename)
#         lengths = {}
#         max_len = 0
#
#         for line in _get_tpx_lines(full_file_path):
#             cols = line.split('\t')
#             if line[0] == '#':
#                 continue
#
#             tfo_start = int(cols[1])
#             tfo_end = int(cols[2])
#             tfo_len = tfo_end - tfo_start
#             lengths[tfo_len] = 1 if tfo_len not in lengths else lengths[tfo_len] + 1
#             if lengths[tfo_len] > max_len:
#                 max_len = lengths[tfo_len]
#
#         # pp.pprint(lengths)
#         plt.xlim(10, 50)
#         plt.ylim(0, max_len + 1000)
#         plt.title("Length distribution")
#         plt.ylabel("#triplexes")
#         plt.xlabel("triplex length")
#
#         sorted_length_keys = sorted(list(lengths.keys())[:50])
#         sorted_length_values = [lengths[x] for x in sorted_length_keys]
#
#         plt.plot(sorted_length_keys, sorted_length_values, label=filename, marker=markers[key], color=colors[key])
#         print lengths
#     plt.legend(loc='upper right')
#     plt.savefig(options.dataOutDir + "/" + out_filename, dpi=150)


def clean_nonmatching_chromosomes(full_file_path):
    """
    Remove triplexes which are located on different chromosomes.
    :param full_file_path:
    :param tpx_content:
    :return:
    """
    print ("Cleaning " + full_file_path)

    tpx_content = _get_tpx_lines(full_file_path)
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


def get_first_line(full_file_path):
    with open(full_file_path) as input_file:
        result = input_file.readlines()[0]
    return result


##########################################################################################################
def create_triplex_superset_tree(full_file_path):
    """

    :param full_file_path:
    :return:
    """
    print ("Creating triplex superset TREE VERSION for file: " + full_file_path)
    tree = intervaltree.IntervalTree()
    content = _get_tpx_lines_by_chromosome(full_file_path)

    output_file = open(full_file_path.replace('.tpx', '.ss.tpx'), 'w')
    output_file.write(get_first_line(full_file_path))

    def compare(iv1, iv2):
        if iv1[0] < iv2[0]:  # start is smaller
            return -1
        elif iv1[0] > iv2[0]:  # start is larger
            return 1
        else:  # starts are equal, check for end point
            return iv2[1] - iv1[1]

    for chrom in content.keys():
        print ("Processing #" + str(len(content[chrom])) + " triplexes from chromosome " + chrom)
        sorted_intervals = set()
        for index, line in enumerate(content[chrom]):
            cols = line.split('\t')
            (tfo_chr, tfo_start_offset) = cols[0].split(':')
            (tts_chr, tts_start_offset) = cols[3].split(':')
            tfo_start_offset = int(tfo_start_offset)
            tts_start_offset = int(tts_start_offset)
            (tfo_start_pos, tfo_end_pos, tts_start_pos, tts_end_pos) = (int(cols[1]) + tfo_start_offset,
                                                                        int(cols[2]) + tfo_start_offset,
                                                                        int(cols[4]) + tts_start_offset,
                                                                        int(cols[5]) + tts_start_offset)
            sorted_intervals.add((tts_start_pos, tts_end_pos, index))
            assert (tfo_chr == tts_chr)

        sorted_intervals = sorted(sorted_intervals, cmp=compare)

        # add first interval always
        tree.add(intervaltree.Interval(sorted_intervals[0][0], sorted_intervals[0][1], sorted_intervals[0][2]))
        output_file.write(content[chrom][sorted_intervals[0][2]])
        for iv in sorted_intervals[1:]:
            overlaps_start = set([x.data for x in tree.top_node.search_overlap([iv[0]])])
            contained = False
            for iv_overlap in tree.top_node.search_overlap([iv[1]]):
                if iv_overlap.data in overlaps_start:
                    contained = True
            if not contained:
                tree.add(intervaltree.Interval(iv[0], iv[1], iv[2]))
                output_file.write(content[chrom][iv[2]])
    output_file.close()


##########################################################################################################
def normalize_sequence_length(full_file_path, nz_len):
    """ Cut each sequence into several sequences of exactly nz_len basepairs. """
    print ("Normalizing sequence length for file: " + full_file_path)
    content = _get_tpx_lines(full_file_path)

    result = []

    for l in range(0, len(content), 2):
        seq_desc = content[l].rstrip()
        seq = content[l + 1].rstrip()
        if len(seq) < nz_len:
            raise ValueError("Found sequence that is smaller than nz_len!: " + seq)
        elif len(seq) > nz_len:
            for shift in range(0, len(seq) - nz_len + 1):
                result.append(seq_desc + "_segment_#" + str(shift) + "\n")
                result.append(seq[shift:shift + nz_len] + "\n")

    output_file = open(full_file_path.replace('.fa', '.nz.fa'), 'w')
    output_file.writelines(result)
    output_file.close()


def merge_overlapping(region_set):
    """

    :param region_set: list of tuples with positions
    :return:
    """
    region_set = sorted(region_set)
    new_region = []
    i = 0
    region_set_len = len(region_set)
    while i < region_set_len - 1:
        start = region_set[i][0]
        end = region_set[i][1]

        j = i
        # j will always be increased by at least 1
        while j < region_set_len and region_set[j][0] <= end:  # start <= end
            end = max(end, region_set[j][1])
            j += 1

        new_region.append((start, end))
        i = j
    return new_region


def jaccard_score(shuffle, ref):
    """

    :param shuffle:
    :param tree: interval tree
    :return:
    """
    score = 0
    # shuffle.merge_overlaps()
    return score


def jaccard(options):
    """
    :param options:
    :return:
    """
    merge_overlapping([])
    return
    if options.jaccardReference is None or options.jaccardQuery is None:
        raise IOError("Jaccard needs ref and query files!")

    ref_region_tree = intervaltree.IntervalTree()
    shuffled_tree = intervaltree.IntervalTree()

    # read in reference
    with open(options.jaccardReference) as ref_file:
        lines = ref_file.readlines()
        for line in lines:
            line = line.split('\t')
            ref_region_tree.add(intervaltree.Interval(int(line[1]), int(line[2])))

    # read query and store
    query_set = []
    with open(options.jaccardQuery) as query_file:
        lines = query_file.readlines()
        for line in lines:
            line = line.split('\t')
            query_set.append((int(line[1]), int(line[2])))
            # shuffled_tree.add(intervaltree.Interval(int(line[1]), int(line[2])))

    # calculate initial jaccard score
    t = time.time()
    print("sort took: " + str(time.time() - t))
    for interval in query_set:
        offset = 0
        shuffled_tree.add(intervaltree.Interval(offset + interval[0], offset + interval[1]))


    init_score = jaccard_score(shuffled_tree, ref_region_tree)


def main(options):
    lazy_import()
    jaccard(options)
    return
    for filename in os.listdir(options.dataInDir):
        if filename.startswith(options.dataPrefix) and filename.endswith(options.dataExtension):
            full_file_path = options.dataInDir + '/' + filename
            if os.path.isfile(full_file_path):
                # clean_nonmatching_chromosomes(full_file_path, tpx_content)
                if options.dataAnalysis == "superset":
                    create_triplex_superset_tree(full_file_path)
                elif options.dataAnalysis == "clean":
                    clean_nonmatching_chromosomes(full_file_path)
                elif options.dataAnalysis == "plot-ld":
                    # plot_length_distributions([full_file_path.replace(".ss.tpx", ".tpx"), full_file_path], [filename.replace(".ss.tpx", ".tpx"), filename], options)
                    plot_length_distributions([full_file_path], [filename], options)
                elif options.dataAnalysis == "normalize":
                    normalize_sequence_length(full_file_path, 20)





