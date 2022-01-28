### input arguments
file_name = "C:/Users/Andrew/Documents/GitHub/bisulfite/example/dsid_C_1461_0.csv"


import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import logsumexp


SCALE = 25


class Dsid:
    '''
    This class is a bisulfite sequence dataset
    '''
    def __init__(self, filename, scale = SCALE):
        '''
        filename: the bisulfite sequence dataset .csv file name
        scale: scale makes all datasets have the same scaled totFeq, thus same variance, 25 is good in practice.
        '''
        self.filename = filename
        self.reads = pd.read_csv(self.filename).set_index("id")
        self.index = [int(x[2:]) for x in self.reads.columns[:-1]]
        self.n = len(self.index)
        self.totSeq = np.zeros(self.n)
        self.totFeq = 0
        for _, row in self.reads.iterrows():
            self.totFeq += row[-1]
            self.totSeq += row[:-1] * row[-1]
        self.totSeq = self.totSeq / self.totFeq * scale
        self.totFeq = scale
        self.a0 = self.totFeq
        self.mu0 = np.sum(self.totSeq) / (self.n * self.totFeq)
        self.alp0 = self.mu0 * self.a0
        self.lam0 = 0
        self.p = []
        self.Yij = np.zeros((self.n, self.n))
        for i in range(self.n):
            self.Yij[i, i:] = self.alp0 + np.cumsum(self.totSeq[i:])
        # filter and backward filter
        self.Pit = np.zeros((self.n, self.n))
        self.Pits = np.zeros((self.n, self.n))
        self.Qjt = np.zeros((self.n, self.n))
        self.Qjts = np.zeros((self.n, self.n))
        self.Bijts = np.zeros((self.n, self.n, self.n))
        self.Bs = np.zeros(self.n)
        self.I = np.zeros(self.n)
        self.Thetat = np.zeros(self.n)
        return

    # beta-binomial
    def pipi(self, a, b, c, d):
        if a == -1 and b == -1:
            return math.lgamma(self.Yij[c][d]) + math.lgamma(self.a0 + (d-c+1)*self.totFeq - self.Yij[c][d]) + math.lgamma(self.a0) - math.lgamma(self.a0 + (d-c+1)*self.totFeq) - math.lgamma(self.alp0) - math.lgamma(self.a0 - self.alp0)
        else:
            return math.lgamma(self.Yij[c][d]) + math.lgamma(self.a0 + (d-c+1)*self.totFeq - self.Yij[c][d]) + math.lgamma(self.a0 + (b-a+1)*self.totFeq) - math.lgamma(self.a0 + (d-c+1)*self.totFeq) - math.lgamma(self.Yij[a][b]) - math.lgamma(self.a0 + (b-a+1)*self.totFeq - self.Yij[a][b])

    def pipipipi(self, a, b, c, d, e, f):
        return self.pipi(a, b, e, f) - self.pipi(-1, -1, c, d)

    def calculatePQ(self):
        for t in range(self.n):
            self.Pits[t, t] = np.log(self.p[t]) + self.pipi(-1, -1, t, t)
            for i in range(t - 1, -1, -1):
                self.Pits[t, i] = np.log(1 - self.p[t]) + self.Pit[t - 1, i] + self.pipi(i, t - 1, i, t)
            lse = logsumexp(self.Pits[t, 0 : t + 1])
            for i in range(t, -1, -1):
                self.Pit[t, i] = self.Pits[t, i] - lse

        for t in range(self.n - 2, -1, -1):
            self.Qjts[t + 1, t + 1] = np.log(self.p[t + 2]) + self.pipi(-1, -1, t + 1, t + 1)
            for j in range(t + 2, self.n):
                self.Qjts[t + 1, j] = self.Qjt[t + 2, j] + self.pipi(t + 2, j, t + 1, j)
            lse = logsumexp(self.Qjts[t + 1, t + 1 :])
            for j in range(t + 1, self.n):
                self.Qjt[t + 1, j] = np.log(1 - self.p[t + 1]) + self.Qjts[t + 1, j] - lse
        return

    def calculateBeta(self):
        for t in range(self.n):
            for i in range(t + 1):
                j = t
                self.Bijts[t, i, j] = np.log(self.p[t + 1]) + self.Pit[t, i]
                for j in range(t + 1, self.n):
                    self.Bijts[t, i, j] = self.Pit[t, i] + self.Qjt[t + 1, j] + self.pipipipi(i, t, t + 1, j, i, j)
            self.Bs[t] = logsumexp(self.Bijts[t, 0 : t + 1, t : self.n].flatten())
        return

    def calculateIt(self):
        self.I[0] = 1
        for t in range(1, self.n - 1):
            self.I[t + 1] = self.p[t + 1] / np.exp(self.Bs[t]) if self.Bs[t] < 500 else 0
        return

    def calculateTheta(self):
        for t in range(self.n):
            for i in range(t + 1):
                for j in range(t, self.n):
                    self.Thetat[t] += np.exp(self.Bijts[t, i, j] - self.Bs[t]) * self.Yij[i, j] / (self.a0 + (j - i + 1) * self.totFeq)
        return

    def choose(self):
        candidate = (0.0, -np.inf)
        lams = -np.log(np.linspace(0.95, 0.9975, 20))
        for lam in lams:
            self.lam0 = lam
            self.p = [1] + [1 - np.e ** (-self.lam0 * (self.index[x] - self.index[x - 1])) for x in range(1, len(self.index))] + [1]
            self.calculatePQ()
            log_likelihood = 0
            for t in range(self.n):
                log_likelihood += logsumexp(self.Pits[t, 0 : t + 1])
            if log_likelihood > candidate[1]:
                candidate = (lam, log_likelihood)
            self.lam0 = 0
            self.p = []
            self.Pit = np.zeros((self.n, self.n))
            self.Pits = np.zeros((self.n, self.n))
            self.Qjt = np.zeros((self.n, self.n))
            self.Qjts = np.zeros((self.n, self.n))

        self.lam0 = candidate[0]
        self.p = [1] + [1 - np.e ** (-self.lam0 * (self.index[x] - self.index[x - 1])) for x in range(1, len(self.index))] + [1]
        return

    def write(self):
        res = pd.DataFrame(np.append(self.I, 0).reshape((1, len(self.index)+1)), columns = self.reads.columns, index = ["prob of change"])
        res1 = pd.DataFrame(np.append(self.Thetat, 0).reshape((1, len(self.index)+1)), columns = self.reads.columns, index = ["prob of response"])
        self.reads = pd.concat([res, res1, self.reads])
        self.reads.to_csv(self.filename[:-4] + "_result.csv")
        return


if __name__ == "__main__":
    sample = Dsid(file_name)
    sample.choose()
    sample.calculatePQ()
    sample.calculateBeta()
    sample.calculateIt()
    sample.calculateTheta()
    sample.write()

    # output posterior estimations
    fig = plt.figure(figsize = (32, 6))
    a = sample.totSeq / sample.totFeq
    b = np.array([math.sqrt(aa * (1 - aa)) / 2 for aa in a])
    p1 = plt.plot(sample.index, sample.Thetat, 'ro', markersize = 10)
    p2 = plt.errorbar(sample.index, a, b, label = 'both limits (default)', fmt = 'bo')
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0", "0.2", "0.4", "0.6", "0.8", "1.0"])
    plt.xlabel("bp locations", fontsize = 20)
    plt.ylabel("frequencies of mutation", fontsize = 20)
    plt.show()

    # output patch size
    t = 0
    peak = sample.totSeq.max()
    while sample.totSeq[t] != peak:
        t += 1
    bij = sample.Bijts[t]
    bij[0 : t + 1, t : sample.n] = np.exp(bij[0 : t + 1, t : sample.n]) / np.exp(sample.Bs[t])
    patch_size, prob = [], []
    for i in range(len(bij)):
        for j in range(len(bij[0])):
            if bij[i, j] > 0.05:
                prob.append(bij[i, j])
                patch_size.append(sample.index[j] - sample.index[i] + 1)
    mean = 0
    aa = sum(prob)
    for i in range(len(prob)):
        prob[i] /= aa
    for i, j in zip(patch_size, prob):
        mean += i * j
    fig = plt.figure(figsize = (10, 10))
    plt.bar(patch_size, prob, width = 0.5)
    plt.xticks(list(range(3, 23)), fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.xlabel("patch size", fontsize = 16)
    plt.ylabel("probability", fontsize = 16)
    fig.suptitle('average size: ' + "{:.2f}".format(mean), fontsize = 24, y = 0.92)
    plt.show()
