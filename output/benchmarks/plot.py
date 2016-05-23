#!/usr/bin/python
import numpy as np
import os
import re
from matplotlib import pyplot as plt
import pprint 

################
# some path defs
PATH_MYERS = "myers/"
PATH_BRUTE = "brute/"

pp = pprint.PrettyPrinter(indent=2)

total_search    = "total_search"

class Value:
    def __init__(self, name):
        self.name   = name
        self.myers    = 0.0
        self.brute    = 0.0

    def __str__(self):
        return str(self.brute) + "," + str(self.myers)


data = {}
calldiff = 0
timediff = 0
cnt = 0
ttime_cs = 0
ttime_ts = 0
ttp_cs_ts= 0

for f in os.listdir(os.getcwd() + "/" + PATH_MYERS):
    if f.endswith(".log"):
        myers_filepath = os.getcwd() + "/" + PATH_MYERS + '/' + f
        brute_filepath = myers_filepath.replace("myers","brute")
        if os.path.isfile(brute_filepath):

            BRUTE_FILE = (open(brute_filepath).read())
            MYERS_FILE = (open(myers_filepath).read())

            #print f
            ts_obj = Value(total_search)
            ts_obj.brute = float((re.search(r'Finished\ program\ within.* ([0-9]*\.[0-9]+)\ seconds', BRUTE_FILE)).group(1))
            ts_obj.myers = float((re.search(r'Finished\ program\ within.* ([0-9]*\.[0-9]+)\ seconds', MYERS_FILE)).group(1))
            #print ts_obj
            #ts_obj.myers = float((re.search(r'triplex\ search\ only.* ([0-9]*\.[0-9]+)\ seconds', MYERS_FILE)).group(1))

            c = int((re.search(r'_c([0-9]+)_', f)).group(1))
            e = int((re.search(r'_e([0-9]+)\.', f)).group(1))
            l = int((re.search(r'_l([0-9]+)_', f)).group(1))

            if l not in data:
                data[l] = {}

            if e not in data[l]:
                data[l][e] = {}

            if c not in data[l][e]:
                data[l][e][c] = {}

            data[l][e][c][total_search] = ts_obj

            cnt += 1

#pp.pprint(data)
###############################################
linestyles = { 1 : '--', 2 : '-', 3 : ':'} 
if True:
    plt.figure(figsize=(10, 10))

    cnt = 211
    for l, l_data in data.iteritems():
        ax = plt.subplot(cnt)
        ax.set_yscale('log')
        plt.ylim(50,25000)

        plt.title("Runtime: brute-force vs. bit-parallel (l : " + str(l) + ')')
        cnt += 1
        plt.xlim(1,25)
        plt.ylabel("Total runtime (s)")
        plt.xlabel("Error rate (%)")
        print "Doing length:  ", l
        errors  = sorted( [v[1] for v in enumerate(l_data) ])

        offset = -2
        for c in range(1,4):
            c_okay = True
            for e, e_data in l_data.iteritems():
                if c not in e_data:
                    c_okay = False
                    break
            if c_okay is False:
                continue

            ts_brute = []
            ts_myers = []
            #for e, e_data in l_data.iteritems():
            for e, e_val in enumerate(errors):
                ts_brute.append(l_data[e_val][c][total_search].brute)
                ts_myers.append(l_data[e_val][c][total_search].myers)



            plt.plot(errors, ts_brute, label="c="+str(c), linestyle=linestyles[c], marker="o", color="red")
            plt.plot(errors, ts_myers, label="c="+str(c), linestyle=linestyles[c], marker="o", color="green")
            for (x,y) in zip(errors, ts_brute):
                    # Annotate the points 5 _points_ above and to the left of the vertex
                ax.annotate('{}'.format(int(y)), xy=(x,y), xytext=(-5, offset), ha='right', textcoords='offset points', color="red")
            
            for (x,y) in zip(errors, ts_myers):
                    # Annotate the points 5 _points_ above and to the left of the vertex
                ax.annotate('{}'.format(int(y)), xy=(x,y), xytext=(+5, offset), ha='left', textcoords='offset points', color="green")

            offset += 5

    plt.legend(loc='upper right')
    plt.savefig("plots/total_search.png", dpi = 100)
    exit()

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

