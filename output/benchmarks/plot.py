#!/usr/bin/python
import numpy as np
import os
import re
from matplotlib import pyplot as plt
import pprint 

################
# some path defs
PATH_MYERS = "./"
PATH_BRUTE = "../mm9/original/"

pp = pprint.PrettyPrinter(indent=2)

total_search    = "total_search"

class Value:
    def __init__(self, name):
        self.name   = name
        self.myers    = 0.0
        self.brute    = 0.0

ts_obj = Value(total_search)

data = {}
calldiff = 0
timediff = 0
cnt = 0
ttime_cs = 0
ttime_ts = 0
ttp_cs_ts= 0

for f in os.listdir(os.getcwd()+"/inverted"):
    if f.endswith(".log") and f.startswith("chr1_serial"):
        if os.path.isfile(os.getcwd()+"/original/"+f):
            ORG_FILE = (open(os.getcwd()+"/original/" + f)).read()
            INV_FILE = (open(os.getcwd()+"/inverted/" + f)).read()

            if re.search(r'seedMultiSeq.* ([0-9]*\.[0-9]+)\ seconds', ORG_FILE) is None:
                continue
            #print f
            ts_obj.brute = float((re.search(r'triplex\ search\ only.* ([0-9]*\.[0-9]+)\ seconds', ORG_FILE)).group(1))
            ts_obj.myers = float((re.search(r'triplex\ search\ only.* ([0-9]*\.[0-9]+)\ seconds', INV_FILE)).group(1))

            io_obj.brute = (re.search(r'IO\ reading.* ([0-9]*\.[0-9]+)\ seconds', ORG_FILE)).group(1)
            io_obj.myers = (re.search(r'IO\ reading.* ([0-9]*\.[0-9]+)\ seconds', INV_FILE)).group(1)

            qf_obj.brute = (re.search(r'qgram-Finder.* ([0-9]*\.[0-9]+)\ seconds', ORG_FILE)).group(1)
            qf_obj.myers = (re.search(r'qgram-Finder.* ([0-9]*\.[0-9]+)\ seconds', INV_FILE)).group(1)

            cs_obj.brute = float((re.search(r'collectSeeds.* ([0-9]*\.[0-9]+)\ seconds', ORG_FILE)).group(1))
            cs_obj.myers = float((re.search(r'collectSeeds.* ([0-9]*\.[0-9]+)\ seconds', INV_FILE)).group(1))

            sm_obj.brute = float((re.search(r'seedMultiSeq.* ([0-9]*\.[0-9]+)\ seconds', ORG_FILE)).group(1))
            sm_obj.myers = float((re.search(r'seedMultiSeq.* ([0-9]*\.[0-9]+)\ seconds', INV_FILE)).group(1))
            
            smc_obj.brute = int((re.search(r'seedMultiSeq was called.* ([0-9]+).*', ORG_FILE)).group(1))
            smc_obj.myers = int((re.search(r'seedMultiSeq was called.* ([0-9]+).*', INV_FILE)).group(1))

            q_threshold = int((re.search(r'_q([0-9]+)_', f)).group(1))
            q = (re.search(r'weight.* ([0-9]+)*', INV_FILE)).group(1)
            c = int((re.search(r'_c([0-9]+)_', f)).group(1))
            e = int((re.search(r'_e([0-9]+)\.', f)).group(1))
            l = int((re.search(r'_l([0-9]+)_', f)).group(1))

            if q not in data:
                data[q] = {}

            if c not in data[q]:
                data[q][c] = {}

            if e not in data[q][c]:
                data[q][c][e] = {}

            data[q][c][e][total_search]     = ts_obj
            data[q][c][e][io_processing]    = io_obj
            data[q][c][e][collect_seeds]    = cs_obj
            data[q][c][e][seed_multi_qgram] = sm_obj

            ttp_cs_ts += (cs_obj.myers / ts_obj.myers) * 100

            calldiff += smc_obj.brute - smc_obj.myers
            timediff += sm_obj.myers - sm_obj.brute
            cnt += 1
# plot settings
print ttp_cs_ts / cnt

print calldiff/cnt
print timediff/cnt
#pp.pprint(data)
###############################################
# seed percentage
if False:
    plt.xlim(1,12)
    plt.ylim(0,1000)
    plt.title("Actual q-gram search function using QGramIndex SA")
    plt.ylabel("# function calls")
    plt.xlabel("q-gram length used for exact matching")

    q_list  = sorted([q[1] for q in enumerate(data) ])
    ts_org  = [data[q][1][5][seed_multi_qgram].brute for q in q_ts]
    ts_inv  = [data[q][1][5][seed_multi_qgram].myers for q in q_ts]

    plt.plot(q_list, ts_org, label="c=1, e=5", linestyle="dashed", marker="o", color="green")
    plt.plot(q_list, ts_inv, label="c=1, e=5", linestyle="dashed", marker="o", color="red")

# ---------------------------------------------
    ts_org  = [data[q][1][10]['ts_o'] for q in q_ts if q < 5]
    ts_inv  = [data[q][1][10]['ts_i'] for q in q_ts if q < 5]
    ts_org += [ts_org[-1]]
    ts_inv += [ts_inv[-1]]


    plt.plot(q_ts, ts_org, label="c=1, e=10", marker="o", color="green")
    plt.plot(q_ts, ts_inv, label="c=1, e=10", marker="o", color="red")
    plt.legend(loc='upper right')
    plt.savefig("plots/total_search.png")


###############################################
# total search
if False:
    plt.xlim(1,12)
    plt.ylim(0,1000)
    plt.title("Total search time without IO")
    plt.ylabel("Total search time without IO (sec)")
    plt.xlabel("Qgram threshold")

    q_ts    = sorted([q[1] for q in enumerate(data) ])
    ts_org  = [data[q][1][5]['ts_o'] for q in q_ts]
    ts_inv  = [data[q][1][5]['ts_i'] for q in q_ts]

    plt.plot(q_ts, ts_org, label="c=1, e=5", linestyle="dashed", marker="o", color="green")
    plt.plot(q_ts, ts_inv, label="c=1, e=5", linestyle="dashed", marker="o", color="red")

# ---------------------------------------------
    ts_org  = [data[q][1][10]['ts_o'] for q in q_ts if q < 5]
    ts_inv  = [data[q][1][10]['ts_i'] for q in q_ts if q < 5]
    ts_org += [ts_org[-1]]
    ts_inv += [ts_inv[-1]]


    plt.plot(q_ts, ts_org, label="c=1, e=10", marker="o", color="green")
    plt.plot(q_ts, ts_inv, label="c=1, e=10", marker="o", color="red")
    plt.legend(loc='upper right')
    plt.savefig("plots/total_search.png")

###############################################
# io
if False:
    plt.clf()
    plt.xlim(0,10)
    plt.ylim(0,150)
    plt.title("IO processing time, without search")
    plt.ylabel("Total IO + processing (sec)")
    plt.xlabel("Maximal # consecutive errors")

    c_io    = [1,2,3]
    io_org  = [data[3][q][5]['io_o'] for q in c_io]
    io_inv  = [data[3][q][5]['io_i'] for q in c_io]

    plt.plot(c_io, io_org, label="q=3, e=5", linestyle="dashed", marker="o", color="green")
    plt.plot(c_io, io_inv, label="q=3, e=5", linestyle="dashed", marker="o", color="red")

# ---------------------------------------------
    c_io    = [1,2,3]
    io_org  = [data[4][q][5]['io_o'] for q in c_io]
    io_inv  = [data[4][q][5]['io_i'] for q in c_io]


    plt.plot(c_io, io_org, label="q=4, e=5", marker="o", color="green")
    plt.plot(c_io, io_inv, label="q=4, e=5", marker="o", color="red")

    plt.savefig("plots/io.png")
