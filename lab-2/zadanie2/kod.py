import numpy as np
import matplotlib.pyplot as plt

def wykres(x, y, nazwa, skala=False):
    if skala:
        plt.yscale('log')
    plt.plot(x, y)
    plt.title(f"{nazwa}")
    #plt.savefig(nazwa)
    plt.show()

def DFT(sygnal):

    N = len(sygnal)
    wynik_tab = np.zeros(N, dtype=complex)

    for k in range(N):
        for n in range(N):
            wynik_tab[k] += sygnal[n] * np.exp((-1j * 2 * np.pi * k * n) / N)

    return wynik_tab

def widmo_amplitudowe(sygnal, dft=True, wbudowane=False):

    if dft == False:
        sygnal_transformata = FFT(sygnal)

    elif dft:
        sygnal_transformata = DFT(sygnal)

    elif wbudowane:
        sygnal_transformata = np.fft.fft(sygnal)

    return np.sqrt(np.real(sygnal_transformata) ** 2 + np.imag(sygnal_transformata) ** 2)[:len(sygnal_transformata) // 2]*2

def skala_decybelowa(sygnal):

    decybele=[]
    for i in sygnal:
        decybele.append(10*np.log(i))

    return decybele

def wykres_x(fs, fi, f, t):
#Funkcja nr 1

    N = int(np.round(t * fs))
    czas = np.linspace(0, t, N)

    def x(t):
       return np.cos(2 * np.pi * f * czas + fi) * np.cos(2.5 * czas ** (0.2) * np.pi)

    wartosci = x(czas)
    #wykres(czas, wartosci, "x")
    return wartosci

def zadanie_2():

    f = 100
    t = 1
    fs = 1000

    czas = np.linspace(0, t, t * fs)
    y = np.sin(2 * np.pi * f * czas)

    widmo = widmo_amplitudowe(y)
    czestotliwosci = np.linspace(0, len(widmo), (t * fs)//2)

    wykres(czestotliwosci, widmo, 'Widmo amplitudowe', skala=False)

    decybele = skala_decybelowa(widmo)
    wykres(czestotliwosci, decybele, 'Skala decybelowa', skala=False)

zadanie_2()

