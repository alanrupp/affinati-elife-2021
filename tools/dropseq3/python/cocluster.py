#!/usr/bin/python

import subprocess
import pandas as pd
import numpy as np

def read_clusters(file):
    d = {}
    with open(file, "r") as f:
        for line in f:
            line.strip()
            if line.strip('\"').startswith("cell"): continue
            else:
                line = line.split(",")
                if line[2] not in d.keys():
                    d[line[2].strip('\"')] = {}
                else:
                    d[line[2].strip('\"')].update({line[0].strip('\"'): line[1].strip('\"')})
    return d

def score(d):
    cells = list(set([k for b in d.keys() for k in d[b].keys()]))
    output = []
    for cell in cells:
        each_score = []
        for other_cell in cells:
            co = []
            for i in list(d.keys()):
                if cell in d[i].keys() and other_cell in d[i].keys():
                    if d[i][cell] == d[i][other_cell]:
                        co.append(1)
                    else:
                        co.append(0)
            if len(co) == 0:
                each_score.append(np.nan)
            else:
                each_score.append(sum(co)/len(co))
        output.append(each_score)
    output = pd.DataFrame(np.array(output))
    output.index = cells
    output.columns = cells
    return output

def rm(file):
    out = f"rm {file}"
    subprocess.call(out, shell=True)
