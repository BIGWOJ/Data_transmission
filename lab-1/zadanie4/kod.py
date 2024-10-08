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

def zadanie_4():
#Funkcja nr 2

    fs = 22050
    t = 1
    #N = int(np.round(t * fs)) ==>> N = fs z powodu, że t = 1
    N = fs
    czas = np.linspace(0, t, N)

    def bk(t, h, k):
        wartosci = (np.sin(np.sin(np.pi * h / 7 * czas) * np.pi * czas * h)) / (2 * h ** 2 + 1)
        return np.transpose(wartosci)

    hk = [2, 5, 25]
    k_tab = [1, 2, 3]

    for numer, h in enumerate(hk):
        for k in k_tab:
            wartosci = bk(t, h, k)
        wykres(czas, wartosci, f"b{numer+1}")

zadanie_4()
