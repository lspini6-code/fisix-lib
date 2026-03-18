import numpy as np
import math

class stat:
    summ = 0.
    summQ = 0.
    N = 0.
    dati = []

    def __init__(self, dati):
        self.dati = dati
        self.summ = sum(self.dati)
        self.summQ = sum([x*x for x in self.dati])
        self.N = len(self.dati)

    def media(self):
        return self.summ/self.N

    def mediana(self):
        return sorted(self.dati)[len(self.dati)//2]

    def var(self, bessel=True):
        vari = self.summQ / self.N - self.media() * self.media()
        if bessel: vari = self.N / (self.N - 1)
        return vari

    def sigma(self, bessel=True):
        return math.sqrt(self.var(bessel))

    def sigma_media(self, bessel=True):
        return math.sqrt(self.var(bessel)/self.N)

    def skew(self):
        media = self.media()
        asimm = 0.
        for x in self.dati:
            asimm += math.pow(x - media, 3)
        asimm = asimm / (self.N * math.pow(self.sigma(), 3))
        return asimm

    def kurt(self):
        media = self.media()
        picc = 0.
        for x in self.dati:
            picc += pow(x - media, 4)
        picc = picc / (self.N * math.pow(self.var(), 2)) - 3.
        return picc

    def append(self, x):
        self.dati.append(x)
        self.summ += x
        self.summQ += x*x
        self.N += 1

    def sturge0(self):
        return int(math.ceil(1. + 3.322*np.log2(self.N)))

    def sturge1(self):
        return int(math.ceil(2.*(self.N**(1/3))))

    def bin_select(self):
        range_ = (min(self.dati), max(self.dati))
        n_bins = self.stuge0()
        bin_conts, bin_ed = np.histogram(self.dati, n_bins, range=range_)
        bin_centres = (bin_ed[:-1] - bin_ed[1:])/2
        return bin_conts, bin_ed, bin_centres
        
