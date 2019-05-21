# proto.py
# Prototyping 2-component EM algorithm considering p-value and LFC.
# Input: tsv file with gene ID, symbol, p-value, and FC.
# Output: ranked list of genes
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture as GM
from sklearn.datasets.samples_generator import make_blobs
#import matplotlib.pyplot as plt


DATA = '5nM_BPA_24H_vs_CNTL_24h.txt'
ID_IND = 1
LFC_IND = -5
P_IND = -1
ids = []
pvals = []
lfcvals = []
allvals = []
with open(DATA, 'r') as infile:
    next(infile)
    for line in infile:
        line = line.strip().split('\t')
        if line[P_IND] == 'NA' or line[LFC_IND] == 'NA':
            pass
        else:
            ids.append(line[ID_IND])
            pvals.append(float(line[P_IND]))
            #lfcvals.append(line[LFC_IND])
            lfcvals.append(abs(float(line[LFC_IND])))
            allvals.append([line[ID_IND], float(line[P_IND]), abs(float(line[LFC_IND]))])
"""
sorted1 = sorted(allvals, key=lambda x: x[1])
with open('sorted_pval.txt', 'w') as outfile:
    for e in sorted1:
        outfile.write('{}\t{}\t{}\n'.format(e[0], e[1], e[2]))
sorted1 = sorted(allvals, key=lambda x: x[2])
with open('sorted_lfc_abs.txt', 'w') as outfile:
    for e in sorted1:
        outfile.write('{}\t{}\t{}\n'.format(e[0], e[1], e[2]))
"""
pvals = np.array(pvals)
lfcvals = np.array(lfcvals)
print(max(lfcvals))
X = []
for i in range(len(ids)):
    X.append([pvals[i], lfcvals[i]])
X = np.array(X)
gmm = GM(n_components=2).fit(X).score_samples(X)
scored = []
for i, score in enumerate(gmm):
    scored.append([ids[i], score, lfcvals[i], pvals[i]])
scored = sorted(scored, key=lambda x: x[1])
with open('results_all_abs.txt', 'w') as outfile:
    outfile.write('gene\tscore\tLFC\tPVAL\n')
    for e in scored:
        outfile.write('{}\t{}\t{}\t{}\n'.format(e[0], e[1], e[2], e[3]))
