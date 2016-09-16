"""

"""

import utils
import os
import re


def compare_result(file1, file2):
    """
    Reads both files via .readlines() in and sorts the lines using `sorted()`.
    :param file1:
    :param file2:
    :return: True if sorted lines match, False otherwise
    """

    with open(file1) as f1:
        with open(file2) as f2:
            c1 = f1.readlines()
            c2 = f2.readlines()
            return sorted(c1) == sorted(c2)


def run_bit_parallel_unit_test(test_id=""):
    """

    :param test_id:
    :return:
    """
    dirs = [test_id] if test_id != "" else [dn for dp, dn, fn in os.walk(utils.PATH_UNIT_TESTS)][0]

    for _dir in dirs:
        # make sure path exists
        if not os.path.isdir(utils.PATH_UNIT_TESTS + _dir):
            print(OSError("Test directory '" + _dir + "' doesn't exist."))
            continue

        test_id = _dir  # get test number / id
        ref_file_name = next((x for x in os.listdir(utils.PATH_UNIT_TESTS + _dir) if x.endswith(".tpx")), None)
        rna_data = utils.PATH_UNIT_TESTS + _dir + "/rna.data"
        dna_data = utils.PATH_UNIT_TESTS + _dir + "/dna.data"

        print ("\n\n##################################################################")
        print ("Running test " + test_id)

        if ref_file_name is None:
            print (IOError("No reference file found in directory for test " + test_id))
            continue

        # get parameters
        L_param = 30
        e_pattern = re.compile("e([0-9]+)")
        e_param = re.search(e_pattern, ref_file_name).group(1)

        c_pattern = re.compile("c([0-9]+)")
        c_param = re.search(c_pattern, ref_file_name).group(1)

        l_pattern = re.compile("l([0-9]+)")
        l_param = re.search(l_pattern, ref_file_name).group(1)

        result_file_name = utils.get_output_filename_from_parameters("bp", l_param, 30, e_param, c_param)
        utils.lib_triplexator.runTriplexator('-ss ' + rna_data + " -ds " + dna_data +
                                             " --bit-parallel -e " + str(e_param) + " -c " + str(c_param) +
                                             " -l " + str(l_param) + " -L " + str(L_param) + " " + " -od " +
                                             utils.PATH_UNIT_TESTS + _dir + " -o " + result_file_name)
        os.remove(utils.PATH_UNIT_TESTS + _dir + "/" + result_file_name + ".log")
        os.remove(utils.PATH_UNIT_TESTS + _dir + "/" + result_file_name + ".summary")

        # evaluate result
        ref_file = utils.PATH_UNIT_TESTS + _dir + "/" + ref_file_name
        result_file = utils.PATH_UNIT_TESTS + _dir + "/" + result_file_name
        print "\n"
        print "Test passed." if compare_result(ref_file, result_file) else "Test failed."
