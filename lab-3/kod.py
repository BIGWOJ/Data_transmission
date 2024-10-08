import numpy as np
import matplotlib.pyplot as plt

#Data rozpoczęcia: 22.04
#Data końcowa: 06.05

def wykres_2(x, y, x2, y2, nazwa, skala=False):
    if skala:
        plt.yscale('log')
    plt.plot(x, y, label="Sygnał zmodulowany")
    plt.plot(x2, y2, label="Ton prosty")
    plt.title(f"{nazwa}")
    plt.legend(loc='upper right')
    plt.savefig(nazwa)
    plt.show()

def wykres(x, y, nazwa, skala=False):
    if skala:
        plt.yscale('log')
    plt.plot(x, y)
    plt.title(f"{nazwa}")
    plt.savefig(nazwa)
    plt.show()

def szerokosc_pasma(decybele_widmo, decybele, czestotliwosci, wykres_nazwa):

    poziom_odnesienia = max(decybele_widmo) - decybele
    obszar_szerokosci = np.where(decybele_widmo >= poziom_odnesienia)
    f_min = np.min(czestotliwosci[obszar_szerokosci])
    f_max = np.max(czestotliwosci[obszar_szerokosci])
    print(f"Wykres: {wykres_nazwa}\nSzerokość pasma {(decybele)}dB: {f_max-f_min}\n")

def ton_prosty(fm, t):
    return np.sin(2*np.pi* fm *t)

def modulacja_amplitudowa(sygnal, ka, fn, t):
    return (ka * sygnal + 1) * np.cos(2*np.pi * fn * t)

def modulacja_fazy(sygnal, kp, f, t):
    return np.cos(2*np.pi * f * t + kp * sygnal)

def modulacja_czestotliwosci(sygnal, kf, fn, fm, t):
    return np.cos(2*np.pi* fn * t + (kf/fm) * sygnal)

def DFT(sygnal):

    N = len(sygnal)
    wynik_tab = np.zeros(N, dtype=complex)

    for k in range(N):
        for n in range(N):
            wynik_tab[k] += sygnal[n] * np.exp((-1j * 2 * np.pi * k * n) / N)

    return wynik_tab

def widmo_amplitudowe(sygnal):
    sygnal_transformata = DFT(sygnal)
    return np.sqrt(np.real(sygnal_transformata) ** 2 + np.imag(sygnal_transformata) ** 2)[:len(sygnal_transformata) // 2]*2

def skala_decybelowa(sygnal):
    return 10*np.log(sygnal)

def zadanie_1_2_3_4():

    #częstotliwość nośna
    fn = 50
    #częsotliwość informacyjna
    fm = 1
    #częstotliwość próbkowania
    fs = 500
    t = 1
    #ilość próbek
    N = fs * t
    czas = np.arange(N) / fs
    sygnal_inf = ton_prosty(fm, czas)
    litery = list('abc')
    decybele = [3,6,12]

#Modulacja amplitudy

    ka_tab = [0.5, 7, 25]

    for ka, litera, decybel in zip(ka_tab, litery, decybele):
        zmodulowany = modulacja_amplitudowa(sygnal_inf, ka, fn, czas)
        wykres_2(x=czas, y=zmodulowany, x2=czas, y2=sygnal_inf, nazwa=f"za_{litera}")

        widmo = widmo_amplitudowe(zmodulowany)
        decybele_widmo = skala_decybelowa(widmo)
        czestotliwosci = np.arange(len(widmo))
        wykres(czestotliwosci, decybele_widmo, f"za_{litera}_widmo")

        szerokosc_pasma(decybele_widmo, decybel, czestotliwosci, f"za_{litera}_widmo")

 #Modulacja fazy

    kp_tab = [0.5, np.pi/2, 5*np.pi]

    for kp, litera, decybel in zip(kp_tab, litery, decybele):
        zmodulowany = modulacja_fazy(sygnal_inf, kp, fn, czas)
        wykres_2(x=czas, y=zmodulowany, x2=czas, y2=sygnal_inf, nazwa=f"zp_{litera}")

        widmo = widmo_amplitudowe(zmodulowany)
        decybele_widmo = skala_decybelowa(widmo)
        czestotliwosci = np.arange(len(widmo))
        wykres(czestotliwosci, decybele_widmo, f"zp_{litera}_widmo")

        szerokosc_pasma(decybele_widmo, decybel, czestotliwosci, f"zp_{litera}_widmo")

 # #Modulacja częstotliwości

    kf_tab = [0.25, np.pi/4, 10*np.pi]

    for kf, litera, decybel in zip(kf_tab, litery, decybele):
        zmodulowany = modulacja_czestotliwosci(sygnal_inf, kf, fn, fm, czas)
        wykres_2(x=czas, y=zmodulowany, x2=czas, y2=sygnal_inf, nazwa=f"zf_{litera}")

        widmo = widmo_amplitudowe(zmodulowany)
        decybele_widmo = skala_decybelowa(widmo)
        czestotliwosci = np.arange(len(widmo))
        wykres(czestotliwosci, decybele_widmo, f"zf_{litera}_widmo")

        szerokosc_pasma(decybele_widmo, decybel, czestotliwosci, f"zf_{litera}_widmo")


zadanie_1_2_3_4()


