import numpy as np
from matplotlib import pyplot as plt

def wykres(x, y, nazwa):
    plt.figure(figsize=(10,5))
    plt.plot(x, y)
    plt.xlabel("Czas")
    plt.ylabel("Amplituda")
    plt.title(f"{nazwa}(t)")
    plt.savefig(nazwa)
    plt.show()

def zadanie_1(fs, fi, f, t):
#Funkcja nr 1

    N = int(np.round(t * fs))
    czas = np.linspace(0, t, N)

    def x(t):
        return np.cos(2 * np.pi * f * t + fi) * np.cos(2.5 * t ** (0.2) * np.pi)

    wartosci = x(czas)
    wykres(czas, wartosci, "x")

fs = 8000
fi = 10
f = 200
t = 2

zadanie_1(fs=fs, fi=fi, f=f, t=t)