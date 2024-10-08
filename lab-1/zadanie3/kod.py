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

def zadanie_3(FS):
#Funkcja nr 1

    t = 1
    #N = int(np.round(t * fs)) ==>> N = fs z powodu, że t = 1
    N = FS
    okres = np.linspace(0, t, N)

    def u(t):
        wartosci = []
        wartosc = 0

        for czas in okres:
            if 0 < czas < 0.1:
                wartosc = np.sin(6* np.pi * czas) * np.cos(5*np.pi * czas)

            elif 0.1 <= czas < 0.4:
                wartosc = -1.1*czas * np.cos(41 * np.pi * (czas**2))

            elif 0.4 <= czas < 0.72:
                wartosc = czas * np.sin(20*(czas**4))

            elif 0.72 <= czas < 1:
                wartosc = 3.3 * (czas - 0.72) * np.cos(27 * np.pi + 1.3)

            wartosci.append(wartosc)

        return np.array(wartosci)

    wartosci = u(t)
    wykres(okres, wartosci, "u")

fs = 8000
fi = 10
f = 200
t = 2

zadanie_3(FS=fs)