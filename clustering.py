### input arguments
file_name = "C:/Users/Andrew/Documents/GitHub/bisulfite/example/dsid_C_1461.csv"
cluster_num = 7


import pandas as pd
from sklearn.cluster import AgglomerativeClustering


def hamming(seq1, seq2):
    '''
    calculate hamming distance between two sequences
    '''
    return sum(s1 != s2 for s1, s2 in zip(seq1, seq2))


def divide(filename, clusters):
    '''
    divide filename into clusters groups and store them separately
    '''
    reads = pd.read_csv(filename).set_index("id")
    seqs = [list(row[:-1]) for _, row in reads.iterrows()]
    sMatrix = [[hamming(seq1, seq2) for seq2 in seqs] for seq1 in seqs]
    labs = AgglomerativeClustering(
        n_clusters = clusters,
        affinity = "precomputed",
        linkage = "complete",
        ).fit(sMatrix).labels_

    for i in range(clusters):
        reads.loc[labs == i].to_csv(filename[:-4] + '_' + str(i) + ".csv")
    return


if __name__ == "__main__":
    divide(file_name, cluster_num)
