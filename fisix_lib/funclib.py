import math
import numpy as np
from rnglib import RandList

'''Classe per la generazione di numeri per la tecnica CrudeMC'''
class AddRec:
    
    def __init__(self, alfa=(np.sqrt(5) - 1)/2):
        self.alfa = alfa
        self.s0 = 0.5
        self.sn = 0.5
    
    def GetNum(self):
        self.sn = (self.sn + self.alfa) % 1
        return self.sn
    
    def SetSeed(self, seed):
        self.s0 = seed
        self.sn = seed
    
    def GetNums(self, n):
        li = []
        for _ in range(n):
            li.append(self.GetNum())
        return li

'''Calcolo dello zero di f in (a, b) con metodo della bisezione'''
def Bisezione(g, a, b, prec=0.0001, it=1000):
    g_min = g(a)
    g_max = g(b)
    if g_min * g_max > 0.:
        raise ValueError("La funzione non può avere estremi di segno uguale")
    for _ in range(it):
        avex = (a + b)/2
        g_ave = g(avex)
        if abs(g_ave) < prec or (b - a)/2 < prec:
            return avex
        if g_ave * g_min > 0.:
            a = avex
            g_min = g_ave
        else:
            b = avex
            g_max = g_ave
    return (a + b)/2.

'''Calcolo dello zero di f in (a, b) con metodo della bisezione ricorsiva'''
def BiRicorsiva(g, a, b, prec=0.0001):
    avex = (a + b)/2.
    if (b - a) < prec:
        return avex
    if g(avex) * g(a) > 0.:
        return BiRicorsiva(g, avex, b, prec)
    else:
        return BiRicorsiva(g, a, avex, prec)

'''Fattoriale di n'''
def Fattoriale(n):
    if (n == 0):
        return 1
    else:
        return n*Fattoriale(n - 1)

'''Calcolo del minimo di f in (a, b) col metodo della sezione aurea'''
def SezAurea(f, a, b, fall=0.0001):
	r = (np.sqrt(5) - 1)/2
	c = b - r*(b - a)
	d = a + r*(b - a)
	fc = f(c)
	fd = f(d)
	while abs(b - a) > fall:
		if fc < fd:
			b, d, fd = d, c, fc
			c = b - r*(b - a)
			fc = f(c)
		else:
			a, c, fc = c, d, fd
			d = a + r*(b - a)
			fd = f(d)	
	return (b + a)/2

'''Tecnica di integrazione Crude Monte Carlo'''
def CrudeMC(g, a, b, punti):
    #Questa tecnica si ricava approssimando E[g(x)]
    rand = RandList(a, b, punti)
    g_rand = [g(x) for x in rand]
    media = np.mean(g_rand)
    var = np.var(g_rand)
    l = b - a
    
    integrale = media * l
    inc_int = np.sqrt(var/punti)*l
    return integrale, inc_int

'''Tecnica di integrazoion Crude Monte Carlo con AddRec (solo tra 0 e 1)'''
def AltCrudeMC(g, punti, alfa):
    seq = AddRec(alfa)
    dense = []
    for _ in range(punti):
        dense.append(seq.GetNum())
    g_dense = [g(x) for x in dense]
    media = np.mean(g_dense)
    media = np.mean(g_dense)
    var = np.var(g_dense)
    l = 1.
    
    integrale = media * l
    inc_int = np.sqrt(var/punti) * l
    return integrale, inc_int

'''Tecnica di integrazione Hit or Miss'''
def HitrMiss(f, a, b, ymax, N):
    x_vals = RandList(a, b, N)
    y_vals = RandList(0., ymax, N)
    
    colpi = 0
    for x, y in zip(x_vals, y_vals):
        if y < f(x): colpi += 1
    
    A = (b - a)*ymax
    p = float(colpi)/float(N) #Frazione di punti che hanno colpito l'area
    I = A*p #Integrale, prodotto dell'area del rettangolo con p
    var = (A**2 * (1 - p) * p) / N #Varianza data dalla Poissoniana (Bernoulli Trial --> Hit OR Miss)
    devst = np.sqrt(var)
    return I, devst, var

'''Tecnica di integrazione Hit or Miss per funzioni 3D'''
def HitrMiss3D(f, xmin, xmax, ymin, ymax, zmax, N):
    colpi = 0
    randX = RandList(xmin, xmax, N)
    randY = RandList(ymin, ymax, N)
    randZ = RandList(0., zmax, N)
    for x, y, z in zip(randX, randY, randZ):
        if z < f(x, y):
            colpi += 1
    Volume = (xmax - xmin)*(ymax - ymin)*zmax
    p = float(colpi)/float(N)
    integrale = Volume * p
    var = (Volume**2 * (1 - p) * p) / N
    devst = np.sqrt(var)
    return integrale, devst, var
