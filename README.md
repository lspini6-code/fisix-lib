# fisix-lib
Libreria Python per calcoli fisici e simulazioni statistiche.

## Moduli
- **rnglib.py**: Generazione di numeri casuali (uniformi, gaussiani, esponenziali, etc.)
- **statlib.py**: Calcoli statistici (media, varianza, deviazione standard, asimmetria, curtosi)
- **funclib.py**: Funzioni matematiche (bisezione, sezione aurea, integrazioni Monte Carlo)
- **luckylib.py**: Funzioni di probabilità (pdf Poisson, esponenziale, Cauchy, likelihood)

## Installazione

### Da GitHub
```bash
pip install git+https://github.com/lspini6-code/fisix-lib.git
```

### Su Google Colab
```python
!pip install git+https://github.com/lspini6-code/fisix-lib.git
from fisix_lib import rnglib, statlib, funclib, luckylib
```

### Da Google Drive (alternativa)
```python
from google.colab import drive
drive.mount('/content/drive')
import sys
sys.path.insert(0, '/content/drive/My Drive/percorso/to/fisix-lib')
from fisix_lib import rnglib, statlib, funclib, luckylib
```

## Utilizzo
Esempio di utilizzo su Colab:
```python
from fisix_lib import rnglib, statlib
import matplotlib.pyplot as plt
# Genera numeri casuali gaussiani
numeri = rnglib.NormBM(mu=0, sigma=1, n=1000)
# Calcola statistiche
stat = statlib.stat(numeri)
print(f"Media: {stat.media()}")
print(f"Deviazione standard: {stat.sigma()}")
# Plot
plt.hist(numeri, bins=30)
plt.show()
```

## Dipendenze
- numpy
- scipy
- matplotlib
- iminuit

## Licenza
MIT
