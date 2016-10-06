import sys
import intervaltree
import random
from matplotlib import pyplot as plt
import pprint

from matplotlib.cbook import is_scalar_or_string

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
    content = utils.get_tpx_lines_by_chromosome(full_file_path)

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
    content = utils.get_tpx_lines(full_file_path)

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
    t = time.time()
    # print("Merging overlapping regions.")
    region_set = sorted(region_set)
    new_region = []
    i = 0
    region_set_len = len(region_set)
    while i < region_set_len:
        start = region_set[i][0]
        end = region_set[i][1]

        j = i
        # j will always be increased by at least 1
        while j < region_set_len and region_set[j][0] <= end:  # start <= end
            end = max(end, region_set[j][1])
            j += 1

        new_region.append((start, end))
        i = j
    # print("merging took: " + str(time.time() - t))
    # print("Region set size before / after merging: " + str(region_set_len) + " / " + str(len(new_region)))
    return new_region


def is_overlapping(t1, t2):
    """

    :param t1: tuple of integers
    :param t2: tuple of integers
    :return: True of t1 and t2 overlap, False otherwise
    """
    return not (t1[0] > t2[1] or t2[0] > t1[1])


def get_overlap(q_int, r_int, q, r):
    """"""
    # no overlap at all
    if q_int[0] > r_int[1]:
        return 0, q, r + 1
    elif q_int[1] < r_int[0]:
        return 0, q + 1, r

    # surely some overlap
    overlap = min(q_int[1], r_int[1]) - max(q_int[0], r_int[0])
    # print ("Overlap: " + str(q_int) + " : " + str(r_int) + " = " + str(overlap))

    if q_int[1] > r_int[1]:
        return overlap, q, r + 1
    else:
        return overlap, q + 1, r


def jaccard_score(query_by_chr, ref_by_chr):
    """

    :param query_by_chr:
    :param ref_by_chr:
    :return:
    """
    intersection = union = 0.0

    for key in ref_by_chr:
        query = merge_overlapping(query_by_chr[key])
        ref = merge_overlapping(ref_by_chr[key])

        query_len = len(query)
        ref_len = len(ref)
        q = r = 0

        # calculate intersection
        while q < query_len and r < ref_len:
            overlap, q, r = get_overlap(query[q], ref[r], q, r)
            intersection += overlap

        # calculate union
        combined = merge_overlapping(query + ref)
        for interval in combined:
            union += interval[1] - interval[0]

    score = intersection / union
    print("Intersection, union, jaccard score: " + str(intersection) + ", " + str(union) + ", " + str(score))
    return score


# TODO cleanup
def jaccard_list(options):
    """
    Perform jaccard test using simple list representations / iterations on the input (no special data structures such
    as interval trees, etc.

    :param options:
    :return:
    """
    print("Running Jaccard test")
    if options.jaccardReference is None or options.jaccardQuery is None:
        raise IOError("Jaccard needs ref and query files!")

    ref_by_chr = {}
    query_by_chr = {}

    # read in reference and store
    with open(options.jaccardReference) as ref_file:
        lines = ref_file.readlines()
        for line in lines:
            line = line.split('\t')
            if line[0] not in ref_by_chr:
                ref_by_chr[line[0]] = []
            ref_by_chr[line[0]].append((int(line[1]), int(line[2])))
    print("Finished reading in reference")

    # read query and store
    with open(options.jaccardQuery) as query_file:
        lines = query_file.readlines()
        for line in lines:
            line = line.split('\t')
            if line[0] not in query_by_chr:
                query_by_chr[line[0]] = []
            query_by_chr[line[0]].append((int(line[1]), int(line[2])))
    print("Finished reading in query")

    # sort input
    t = time.time()
    ref_by_chr = {k: sorted(v) for k, v in ref_by_chr.items()}
    query_by_chr = {k: sorted(v) for k, v in query_by_chr.items()}
    print ("sort took: " + str(time.time() - t))

    # calculate initial jaccard score
    init_score = jaccard_score(query_by_chr, ref_by_chr)
    print("Jaccard initial score: " + str(round(init_score, 4)))

    # dictionary with maximum length for each chromosome
    max_chr_len = {k: max([interval[1] for interval in v]) for k, v in query_by_chr.items()}

    # do n random permutations
    n = 10
    scores = []
    for i in range(n):
        r = random.SystemRandom()
        random_query_by_chr = {
            k: [(random_offset - (interval[1] - interval[0]), random_offset) for (random_offset, interval)
                in ((r.randint(0, max_chr_len[k]), interval) for interval in v)]
            for k, v in query_by_chr.items()
            }
        score = round(jaccard_score(random_query_by_chr, ref_by_chr), 4)
        print("Jaccard random score #" + str(n) + ": " + str(score))
        scores.append(score)

    p = float(len([x for x in scores if x > init_score])) / n
    print ("Joseph's p ([] > real / nr_tests): " + str(p))


def plot_len_dist_multi(options):
    """
    Plot normalized length distribution of triplexes within multiple regions against triplexes outside these regions.
    Input must be .bed files!
    :param options:
    :return:
    """
    utils.check_io_files([options.lenDistTpx, options.lenDistRegions])

    region_files = options.lenDistRegions.split(' ')
    # region_names = [ for f in options.lenDistRegions.split(' ')]

    # ###############################
    # read in triplexes by chromosome
    tpx_chrom_dict = {}  # each key contains an interval tree
    tpx_used = []  # keep track which triplexes have been used
    with open(options.lenDistTpx) as tpx_file:
        index = 0
        for tpx in tpx_file.readlines():
            chrom, start, end = tpx.split('\t')[0], int(tpx.split('\t')[1]), int(tpx.split('\t')[2])
            if chrom not in tpx_chrom_dict:
                tpx_chrom_dict[chrom] = intervaltree.IntervalTree()

            tpx_used.append(False)
            tpx_chrom_dict[chrom].add(intervaltree.Interval(start, end, index))
            index += 1

    # ###############################
    # init dictionary containing length distributions by region
    dist_by_region = {}
    # iterate through each region and populate distributions
    for r in region_files:
        region_name = r.split("/")[-1]
        dist_by_region[region_name] = {}
        # open region file
        with open(r) as region_file:
            for region in region_file.readlines():
                chrom, start, end = region.split('\t')[0], int(region.split('\t')[1]), int(region.split('\t')[2])
                if chrom not in tpx_chrom_dict:
                    continue
                # ###############################
                # check if this region overlaps / contains any triplex
                overlapped_tpxs = tpx_chrom_dict[chrom].top_node.search_overlap([start, end])

                for interval in overlapped_tpxs:
                        # if overlapping and not used, increase count
                        index = interval.data
                        if tpx_used[index] is False:
                            tpx_used[index] = True

                        if interval.length() not in dist_by_region[region_name]:
                            dist_by_region[region_name][interval.length()] = 0
                        dist_by_region[region_name][interval.length()] += 1
                        # print(chrom, tpx, " == ", start, end)
                        tpx_used[index] = True

    # ###############################
    # normalize distributions
    for region_name, dist in dist_by_region.items():
        s = float(sum(dist.values()))
        for l in dist:
            dist[l] /= s

    # ###############################
    # plot distributions
    colors = ["green", "red", "blue", "yellow"]
    key = 0
    for region_name, dist in dist_by_region.items():
        keys = []
        vals = []
        for k in range(10, 100):
            keys.append(k)
            vals.append(dist[k] if k in dist else 0.0)

        plt.xlim(10, 100)
        plt.ylim(0, 1)
        plt.title("Length distribution")
        plt.ylabel("#triplexes")
        plt.xlabel("triplex length")
        plt.plot(keys, vals, label=region_name, color=colors[key])
        key += 1
    plt.show()


def main(options):
    lazy_import()

    if options.dataAnalysis == "len-dist-multi":
        plot_len_dist_multi(options)


    # for filename in os.listdir(options.dataInDir):
    #     if filename.startswith(options.dataPrefix) and filename.endswith(options.dataExtension):
    #         full_file_path = options.dataInDir + '/' + filename
    #         if os.path.isfile(full_file_path):
    #             # clean_nonmatching_chromosomes(full_file_path, tpx_content)
    #             if options.dataAnalysis == "superset":
    #                 create_triplex_superset_tree(full_file_path)
    #             elif options.dataAnalysis == "clean":
    #                 clean_nonmatching_chromosomes(full_file_path)
    #             elif options.dataAnalysis == "plot-ld":
    #                 # plot_length_distributions([full_file_path.replace(".ss.tpx", ".tpx"), full_file_path], [filename.replace(".ss.tpx", ".tpx"), filename], options)
    #                 plot_length_distributions([full_file_path], [filename], options)
    #             elif options.dataAnalysis == "normalize":
    #                 normalize_sequence_length(full_file_path, 20)
