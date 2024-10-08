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
        return np.cos(2 * np.pi * f * czas + fi) * np.cos(2.5 * czas ** (0.2) * np.pi)

def zadanie_2(fs, fi, f, t):
#Zestaw funkcji nr 1

    N = int(np.round(t * fs))
    czas = np.linspace(0, t, N)

    def x(t):
       return np.cos(2 * np.pi * f * czas + fi) * np.cos(2.5 * czas ** (0.2) * np.pi)

    wartosci = x(t)
    wykres(czas, wartosci, "x")

    def y(t):
        return (x(czas)* czas) / (3+ np.cos(20* np.pi * czas))

    wartosci = y(t)
    wykres(czas, wartosci, "y")

    def z(t):
        return czas**2* np.abs(x(t) * y(czas) - (2/ 10+y(czas)))

    wartosci = z(t)
    wykres(czas, wartosci, "z")

    def v(t):
        return z(czas)**3 + 3*np.sin(z(czas) * y(czas)) * np.abs(y(czas) - x(czas))

    wartosci = v(t)
    wykres(czas, wartosci, "v")

fs = 8000
fi = 10
f = 200
t = 2

zadanie_2(fs=fs, fi=fi, f=f, t=t)