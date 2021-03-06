fn = "C:/Users/Yingru/Source/Repos/bisf/data/dsid_C_1461.csv"
cn = 7


import pandas as pd
from sklearn.cluster import AgglomerativeClustering


def hamming(seq1, seq2):
    return sum(s1 != s2 for s1, s2 in zip(seq1, seq2));


def divide(filename, c):
    reads = pd.read_csv(filename).set_index("id");
    seqs = [list(line[1][:-1]) for line in reads.iterrows()];
    sMatrix = [[hamming(seq1, seq2) for seq2 in seqs] for seq1 in seqs];
    labs = AgglomerativeClustering(n_clusters = c, affinity = "precomputed", linkage = "complete").fit(sMatrix).labels_;
    for i in range(c):
        reads.iloc[labs == i].to_csv(filename[:-4] + "_" + str(i) + ".csv");
    return;


if __name__ == "__main__":
    divide(fn, cn);