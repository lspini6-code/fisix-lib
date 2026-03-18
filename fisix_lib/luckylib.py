import numpy as np
from funclib import Fattoriale, Bisezione

'''pdf Poissoniana'''
def PoissPdf(x, Lambda):
    return (Lambda**int(x)) * np.exp(-Lambda) / Fattoriale(conts)

'''pdf Esponenziale'''
def ExpPdf(x, Lambda):
    return np.exp(-x * Lambda) * Lambda

'''pdf Gaussiana asimmetrica'''
def DoubleGauss(x, my, sigmaS, sigmaD):
    A = 2./(np.sqrt(2.*np.pi) * (sigmaS + sigmaD))
    if x < mu:
        return A * np.exp(-0.5 * ((x - mu)/sigmaS)**2)
    else:
        return A * np.exp(-0.5 * ((x - mu)/sigmaD)**2)

'''pdf di Cauchy'''
def CauchyPdf(x, M, gamma):
    return gamma/(np.pi*((x - M)**2 + gamma**2))

'''Likelihood, funzione di verosimiglianza'''
def Like(th, pdf, dati):
    p = pdf(np.arra(dati), th)
    return np.prod(p)

'''Likelihood della pdf esponenziale, poco utile'''
def ExpLike(Lambda, dati):
    n = len(dati)
    return Lambda**n * np.exp(-np.sum(dati)*Lambda)

'''Log-Likelihood, logaritmo della funzione di verosimiglianza'''
def LogLike(th, pdf, dati):
    p = pdf(np.array(dati), th)
    p = np.where(p > 0., p, 1e-300)
    return np.sum(np.log(p))

'''Negativo della Log-LikeLikelihood, utile per trovare il massimo della Likelihood'''
def NegLogLike(th, pdf, dati):
    return -LogLike(th, pdf, dati)

'''Log-LL della pdf esponenziale'''
def ExpLogLike(Lambda, dati):
    n = len(dati)
    return n*np.log(Lambda) - np.sum(dati) * Lambda

'''Log-LL della pdf Poissoniana'''
def PoissLogLike(Lambda, conts):
    return np.sum(np.log(np.array([PoissPdf(int(k), Lambda) for k in conts])))

'''Metodo della sezione aurea per il calcolo del massimo della NegLog-LL'''
def SezAureaLL(g, pdf, dati, a, b, prec=0.0001):
	r = (np.sqrt(5) - 1)/2
	x2 = x1 - r*(x1 - x0)
	x3 = x0 + r*(x1 - x0)
	g2 = g(x2, pdf, eventi)
	g3 = g(x3, pdf, eventi)
	
	while abs(x1 - x0) > prec:
		if g2 < g3:
			x1 = x3
			x3 = x2
			g3 = g2
			x2 = x1 - r*(x1 - x0)
			g2 = g(x2, pdf, eventi)
		else:
			x0 = x2
			x2 = x3
			g2 = g3
			x3 = x0 + r*(x1 - x0)
			g3 = g(x3, pdf, eventi)
	
	return (x0 + x1)/2

'''Stima grafica delle incertezze'''
def GraphicInc(pdf, dati, a, b):
    f = lambda x: LogLike(x, pdf, dati)
    th_hat = SezAurea(-LogLike, pdf, dati, a, b)
    LL_max = f(th_hat)
    Altf = lambda x: f(x) - LL_max + 0.5
    Left = Bisezione(Altf, a, th_hat)
    Right = Bisezione(Altf, th_hat, b)
    return th_hat, Left, Right

'''Funzione per trovare Lambda ± sigma_Lambda'''
def intLL(g, pdf, dati, a, b, ylev, Lamb, prec=0.0001):
    def gpr(x):
        return g(x, pdf, dati, Lamb) - ylev
    avex = a
    while (b - a) > prec:
        avex = (b + a)/2.
        if gpr(avex) * gpr(avex) > 0. : a = avex
        else: b = avex
    return avex

'''Calcola la frazione di elementi di Q2list che sono minori di Q2th'''
def CalcProb(Q2list, Q2th):
    return len([Q2 for Q2 in Q2list if Q2 < Q2th])/len(Q2list)

'''Rapporto tra Log-LL di Lambda e Log-LL di Lambda che massimizza Log-LL'''
def LLR(lambval, pdf, dati, lambmax):
    llval = LogLike(lambval, pdf, dati)
    llmax = LogLike(lambmax, pdf, dati)
    return llval - llmax

#Trascrive i dati in un file in liste.
'''
def LeggiFile():
    with open('dati.txt', 'r') as doc:
        stuff = [line.strip().split() for line in doc]
    cols = list(zip(*stuff))
    c1 = [float(x) for x in cols[0]]
    c2 = [float(x) for x in cols[1]]
    c3 = [float(x) for x in cols[2]]
    return c1, c2, c3
'''
