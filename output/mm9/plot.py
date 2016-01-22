#!/usr/bin/python
import numpy as np
import os
import re
from matplotlib import pyplot as plt

TS_ORG = 0
TS_INV = 0

IO_ORG = 0
IO_INV = 0

QF_ORG = 0
QF_INV = 0

CS_ORG = 0
CS_INV = 0

data = {}

for f in os.listdir(os.getcwd()+"/inverted"):
    if f.endswith(".log"):
        if os.path.isfile(os.getcwd()+"/original/"+f):
            ORG_FILE = (open(os.getcwd()+"/original/" + f)).read()
            INV_FILE = (open(os.getcwd()+"/inverted/" + f)).read()

            TS_ORG = (re.search(r'triplex\ search\ only.* ([0-9]*\.[0-9]+)\ seconds', ORG_FILE)).group(1)
            TS_INV = (re.search(r'triplex\ search\ only.* ([0-9]*\.[0-9]+)\ seconds', INV_FILE)).group(1)

            IO_ORG = (re.search(r'IO\ reading.* ([0-9]*\.[0-9]+)\ seconds', ORG_FILE)).group(1)
            IO_INV = (re.search(r'IO\ reading.* ([0-9]*\.[0-9]+)\ seconds', INV_FILE)).group(1)

            QF_ORG = (re.search(r'qgram-Finder.* ([0-9]*\.[0-9]+)\ seconds', ORG_FILE)).group(1)
            QF_INV = (re.search(r'qgram-Finder.* ([0-9]*\.[0-9]+)\ seconds', INV_FILE)).group(1)

            CS_ORG = (re.search(r'collectSeeds.* ([0-9]*\.[0-9]+)\ seconds', ORG_FILE)).group(1)
            CS_INV = (re.search(r'collectSeeds.* ([0-9]*\.[0-9]+)\ seconds', INV_FILE)).group(1)
            
            q = int((re.search(r'_q([0-9]+)_', f)).group(1))
            c = int((re.search(r'_c([0-9]+)_', f)).group(1))
            e = int((re.search(r'_e([0-9]+)\.', f)).group(1))
            l = int((re.search(r'_l([0-9]+)_', f)).group(1))

            if q not in data:
                data[q] = {}

            if c not in data[q]:
                data[q][c] = {}

            if e not in data[q][c]:
                data[q][c][e] = {}

            data[q][c][e]['ts_o'] = TS_ORG
            data[q][c][e]['ts_i'] = TS_INV

            data[q][c][e]['io_o'] = IO_ORG
            data[q][c][e]['io_i'] = IO_INV

# plot settings

###############################################
# total search
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
# Show result on screen
