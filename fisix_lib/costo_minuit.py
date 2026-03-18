'''
Il programma confronta il valore del Q^2 al variare della funzione di costo
in particolare si considerano BinnedNLL, ExtendedBinnedNLL e LeastSquares.
'''
import numpy as np
import matplotlib.pyplot as plt
from statlib import stat
from rnglib import NormTCL
from math import sqrt, ceil
from iminuit import Minuit
from iminuit.cost import LeastSquares, ExtendedBinnedNLL, BinnedNLL
from scipy.stats import norm, chi2

#Modello per BinnedNLL, cdf della gaussiana normalizzata.
def ModBinLL(bin_ed, mu, sigma):
    return norm.cdf(bin_ed, mu, sigma)

#Modello per ExtendedBinnedNLL, cdf della gaussiana. Tiene in cosiderazione il numero di eventi.
def ModBinELL(bin_ed, n_sig, mu, sigma):
    return n_sig * norm.cdf(bin_ed, mu, sigma)

#Modello per LeastSquares, pdf della gaussiana perché il confronto è tra la pdf e la forma dell'istogramma.
def ModBinLS(x, n_evt, mu, sigma, bin_width):
    return n_evt * norm.pdf(x, mu, sigma) * bin_width

def Sturges(n):
    return int(ceil(1. + 3.322*np.log2(n)))

toys = 5000 #Numero di esperimenti.
n_mis = 500 #Numero di misure per esperimento.
QLL, QELL, QLS = [], [], [] #Liste di valori di Q^2.
ndofLL, ndofELL, ndofLS = 0, 0, 0 #Inizializzazione numero di gradi di libertà.
n_bins = Sturges(n_mis) #Numero di bins per il fit.
N_bins = Sturges(toys) #Numero di bins per l'istogramma.

for _ in range(toys):
    eve = NormTCL(1., 0.4, n_mis) #Genera n_mis numeri casuali secondo una gaussiana con media 1 e dev.std. 0.4.
    my_stat = stat(eve)
    media = my_stat.media() #Media di eve.
    devst = my_stat.sigma() #Dev.std. di eve.
    range_ = (min(eve), max(eve)) #Range in cui costruisco l'istogramma.
    bin_conts, bin_ed = np.histogram(eve, n_bins, range=range_) #bin_conts --> frequenze nei bin, bin_ed --> intervalli bin.
    n_evt = sum(bin_conts) #Somma delle frequenze: numero di eventi totali.
    bin_centres = 0.5 * (bin_ed[1:] + bin_ed[:-1]) #Centri dei bin.
    bin_width = bin_centres[1] - bin_centres[0] #Larghezza dei bin.
    sigma_cont = [sqrt(n) for n in bin_conts] #La dev.std. così definita deriva dalla dev.std. della Poissoniana.
    
    '''Minimizzazione con BinnedNLL'''
    costLL = BinnedNLL(bin_conts, bin_ed, ModBinLL) #Sintassi: BinnedNLL(osservati, intervalli, modello)
    my_binLL = Minuit(costLL, mu=media, sigma=devst) #Inizializzo mu e sigma per facilitare minuit.
    my_binLL.limits['sigma'] = (0, None) #Impongo sigma maggiore di 0.
    my_binLL.migrad() #Minimizzazione.
    my_binLL.hesse() #Calcolo delle incertezze.
    if my_binLL.valid:
        QLL.append(my_binLL.fval)
        ndofLL = my_binLL.ndof
    
    '''Minimizzazione con ExtendedBinnedNLL'''
    costELL = ExtendedBinnedNLL(bin_conts, bin_ed, ModBinELL) #Sintassi molto simile a BinnedNLL.
    my_binELL = Minuit(costELL, n_sig=n_evt, mu=media, sigma=devst) #Inizializzo n_sig con n_evt.
    my_binELL.limits['n_sig', 'sigma'] = (0, None)
    my_binELL.migrad()
    my_binELL.hesse()
    if my_binELL.valid:
        QELL.append(my_binELL.fval)
        ndofELL = my_binELL.ndof
    
    '''Minimizzazione con LeastSquares'''
    costLS = LeastSquares(bin_centres, bin_conts, sigma_cont, ModBinLS) #Sintassi: x-->bin_centres, y-->bin_conts, sigma-->sigma_cont, modello.
    my_binLS = Minuit(costLS, n_evt=n_evt, mu=media, sigma=devst, bin_width=bin_width)
    my_binLS.fixed['n_evt', 'bin_width'] = True #Escludo dal fit quesi due parametri.
    my_binLS.limits['sigma'] = (0, None)
    my_binLS.migrad()
    my_binLS.hesse()
    if my_binLS.valid:
        QLS.append(my_binLS.fval)
        ndofLS = my_binLS.ndof

statLL, statELL, statLS = stat(QLL), stat(QELL), stat(QLS)
print("Media con BinnedNLL: ", statLL.media())
print("Media con ExtendedBinnedNLL: ", statELL.media())
print("Media con LeastSquares: ", statLS.media())

'''Visualizzazione grafica con matplotlib'''
fig, axes = plt.subplots(nrows=3, ncols=1)
axes[0].set_title('Istogramma di Q^2 al variare della funzione di costo.')

axes[0].set(xlabel='Q^2 BinnedNLL', ylabel='Eventi per bin')
axes[0].hist(QLL, N_bins, color='red')

axes[1].set(xlabel='Q^2 ExtendedBinnedNLL', ylabel='Eventi per bin')
axes[1].hist(QELL, N_bins, color='orange')

axes[2].set(xlabel='Q^2 LeastSquares', ylabel='Eventi per bin')
axes[2].hist(QLS, N_bins, color='yellow')

plt.show()
