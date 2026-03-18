import random
import math
import numpy as np
from scipy.stats import norm

'''Signolo numero casuale (uniforme) tra a e b'''
def RandRange(a, b):
    return a + random.random() * (b - a)

'''Numero intero casuale (uniforme) tra a e b'''
def IntRange(a, b):
    return math.ceil(RandRange(a, b))

'''Lista di n numeri casuali (uniformi) tra a e b'''
def RandList(a, b, n):
    r = []
    for _ in range(n):
        r.append(RandRange(a, b))
    return r

'''Lista di n numeri interi casuali (uniformi) tra a e b'''
def RandIntList(a, b, n):
    r = []
    for _ in range(n):
        r.append(IntRange(a, b))
    return r

'''Funzione di Box-Mueller'''
def BoxMue(x1, x2):
    return np.sqrt(-2.*np.log(x1)) * np.cos(2.*np.pi*x2)

'''Generazione di n numeri con la tecnica Try and Catch'''
def TrynCatch(f, a, b, ymax, n):
    r = []
    for _ in range(n):
        x = RandRange(a, b)
        y = RandRange(0., ymax)
        while y > f(x):
            x = RandRange(a, b)
            y = RandRange(0., ymax)
        r.append(x)
    return r

'''Generazione di n numeri secondo il Teorema Centrale del Limite'''
def CentralLimit(a, b, n, somme=500):
    r = []
    for _ in range(n):
        y = 0.
        for _ in range(somme):
            y += RandRange(a, b)
        r.append(y/somme)
    return r

'''Numeri gaussiani generati dalla Box-Mueller'''
def NormBM(mu, sigma, n):
    r = []
    for _ in range(n):
        x1 = RandRange(0.0001, 1.)
        x2 = RandRange(0., 1.)
        r.append(sigma*BoxMue(x1, x2) + mu)
    return r

'''Numeri gaussiani generati col Try and Catch'''
def NormTAC(mu, sigma, n):
    r = []
    ymax = norm.pdf(mu, mu, sigma)
    a = mu - 4*sigma
    b = mu + 4*sigma
    for _ in range(n):
        x = RandRange(a, b)
        y = RandRange(0., ymax)
        while y > norm.pdf(x, mu, sigma):
            x = RandRange(a, b)
            y = RandRange(0., ymax)
        r.append(x)
    return r

'''Numeri gaussiani generati col Teorema Centrale del Limite'''
def NormTCL(mu, sigma, toys, n=500):
    medie = []
    a = mu - (sigma * math.sqrt(12.*n))/2.
    b = mu + (sigma * math.sqrt(12.*n))/2.
    for _ in range(toys):
        x = RandList(a, b, n)
        media = np.mean(x)
        medie.append(media)
    return medie

'''Funzione inversa della cdf esponenziale'''
def InverseExpCdf(y, Lambda):
    return -1./Lambda * np.log(1. - y)

'''Numero casuale che segue una pdf esponenziale'''
def ExpNum(Lambda):
    return -1./Lambda * np.log(1. - random.random())

'''Lista di nueri casuali distribuiti secondo una pdf esponenziale, generati col metodo della fumnzione inversa'''
def ExpList(Lambda, n):
    r = []
    for _ in range(n):
        x = InverseExpCdf(random.random(), Lambda)
        r.append(x)
    return r

'''Generazione di un punto (l, alfa) con l r.v. gaussiana e alfa r.v. uniforme'''
def RandomStep(mu, sigma):
    passi = 1
    lungh = NormTAC(mu, sigma, passi)[0]
    angolo = RandRange(0., 2.*np.pi)
    x = lungh * np.cos(angolo)
    y = lungh * np.sin(angolo)
    return a, y

'''Lista di punti (l, alfa)'''
def RandomWalk(mu, sigma, passi):
    lungh = NormTAC(mu, sigma, passi)
    angoli = RandList(0., 2.*np.pi, passi)
    x_coord, y_coord = [0.], [0.]
    for l, alfa in zip(lungh, angoli):
        x_coord.append(l*np.cos(alfa))
        y_coord.append(l*np.sin(alfa))
    X, Y = sum(x_coord), sum(y_coord)
    return X, Y, x_coord, y_coord

'''Simulazione di eventi Poissoniani. Sfurtta il fatto che i tempi seguono una pdf esponenziale'''
def GenPoiss(Lambda):
    tm = 1.
    conteggi = []
    N = 10000
    for _ in range(N):
        tem, con = 0., 0.
        while tem <= tm:
            t = ExpNum(1./Lambda)
            tem += t
            if tem <= tm:
                con += 1
        conteggi.append(con)
    return conteggi
