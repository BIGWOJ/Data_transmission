import numpy as np
import matplotlib.pyplot as plt
import time

def wykres(x, y, nazwa, skala=False):
    if skala:
        plt.yscale('log')
    plt.plot(x, y)
    plt.title(f"{nazwa}")
    plt.savefig(nazwa)
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

def FFT(sygnal):
    sygnal_transformata = np.fft.fft(sygnal)
    return np.sqrt(np.real(sygnal_transformata) ** 2 + np.imag(sygnal_transformata) ** 2)[:len(sygnal_transformata) // 2]*2

def wykres_x(fs, fi, f, t):
#Funkcja nr 1

    N = int(np.round(t * fs))
    czas = np.linspace(0, t, N)

    def x(t):
       return np.cos(2 * np.pi * f * czas + fi) * np.cos(2.5 * czas ** (0.2) * np.pi)

    wartosci = x(czas)
    #wykres(czas, wartosci, "x")
    return wartosci

def wykres_yzv(fs, fi, f, t):
#Zestaw funkcji nr 1

    N = int(np.round(t * fs))
    czas = np.linspace(0, t, N)

    def x(t):
       return np.cos(2 * np.pi * f * czas + fi) * np.cos(2.5 * czas ** (0.2) * np.pi)

    wartosci_x = x(t)
    #wykres(czas, wartosci, "x")

    def y(t):
        return (x(czas)* czas) / (3+ np.cos(20* np.pi * czas))

    wartosci_y = y(t)
    #wykres(czas, wartosci, "y")

    def z(t):
        return czas**2* np.abs(x(t) * y(czas) - (2/ 10+y(czas)))

    wartosci_z = z(t)
    #wykres(czas, wartosci, "z")

    def v(t):
        return z(czas)**3 + 3*np.sin(z(czas) * y(czas)) * np.abs(y(czas) - x(czas))

    wartosci_v = v(t)
    #wykres(czas, wartosci, "v")

    return wartosci_y, wartosci_z, wartosci_v

def wykres_u(FS):
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
    #wykres(okres, wartosci, "u")
    return wartosci

def wykres_b(wbudowane=False, dft=True, wykresy=True, pomiar_czasu=False):

    fs = 2000
    t = 1
    # N = int(np.round(t * fs)) ==>> N = fs z powodu, że t = 1
    N = fs
    czas = np.linspace(0, t, N)

    def bk(t, h, k):
        wartosci = (np.sin(np.sin(np.pi * h / 7 * czas) * np.pi * czas * h)) / (2 * h ** 2 + 1)
        return np.transpose(wartosci)

    hk = [2, 5, 25]
    k_tab = [1, 2, 3]
    #wartosci_tab = []

    for numer, h in enumerate(hk):
        for k in k_tab:
            wartosci = bk(t, h, k)
        #wartosci_tab.append(wartosci)

        if dft == False:
            czas_poczatek_sygnal = 0
            czas_poczatek_sygnal = time.time()
            widmo = widmo_amplitudowe(wartosci, dft=False)
            if pomiar_czasu:
                print(f"Czas FFT dla sygnału b{numer + 1}: {(time.time() - czas_poczatek_sygnal):.15f}")


        if dft:
            czas_poczatek_sygnal = time.time()
            widmo = widmo_amplitudowe(wartosci)
            if pomiar_czasu:
                print(f"Czas DFT dla sygnału b{numer + 1}: {time.time() - czas_poczatek_sygnal}")

        czestotliwosci = np.arange(len(widmo))

        if wykresy:
            wykres(czestotliwosci, widmo, f"b{numer + 1}", skala=True)

        if wbudowane:
            widmo = widmo_amplitudowe(wartosci)
            wykres(czestotliwosci, widmo, f"b{numer + 1}", skala=True)
            widmo = widmo_amplitudowe(wartosci, wbudowane=True)
            wykres(czestotliwosci, widmo, f"b{numer + 1} - wbudowane fft", skala=True)

    #return wartosci_tab

def zadanie_4():

    fs = 500
    fi = 10
    f = 200
    t = 1

    x = wykres_x(fs, fi, f, t)
    czas_poczatek = time.time()
    czas_poczatek_sygnal = time.time()
    widmo = widmo_amplitudowe(x, dft=False)
    print(f"Czas FFT dla sygnału x: {(time.time() - czas_poczatek_sygnal):.15f}")

    y, z, v = wykres_yzv(fs, fi, f, t)
    czas_poczatek_sygnal = time.time()
    widmo = widmo_amplitudowe(y, dft=False)
    print(f"Czas FFT dla sygnału y: {(time.time() - czas_poczatek_sygnal):.15f}")

    czas_poczatek_sygnal = time.time()
    widmo = widmo_amplitudowe(z, dft=False)
    print(f"Czas FFT dla sygnału z: {(time.time() - czas_poczatek_sygnal):.15f}")

    czas_poczatek_sygnal = time.time()
    widmo = widmo_amplitudowe(v, dft=False)
    print(f"Czas FFT dla sygnału v: {(time.time() - czas_poczatek_sygnal):.15f}")

    u = wykres_u(fs)
    czas_poczatek_sygnal = time.time()
    widmo = widmo_amplitudowe(u, dft=False)
    print(f"Czas FFT dla sygnału u: {(time.time() - czas_poczatek_sygnal):.15f}")

    wykres_b(dft=False, wykresy=False, pomiar_czasu=True)

    czas_fft = time.time() - czas_poczatek


    # ========================


    x = wykres_x(fs, fi, f, t)
    czas_poczatek = time.time()
    czas_poczatek_sygnal = time.time()
    widmo = widmo_amplitudowe(x)
    print(f"\nCzas DFT dla sygnału x: {time.time() - czas_poczatek_sygnal}")

    y, z, v = wykres_yzv(fs, fi, f, t)
    czas_poczatek_sygnal = time.time()
    widmo = widmo_amplitudowe(y)
    print(f"Czas DFT dla sygnału y: {time.time() - czas_poczatek_sygnal}")

    czas_poczatek_sygnal = time.time()
    widmo = widmo_amplitudowe(z)
    print(f"Czas DFT dla sygnału z: {time.time() - czas_poczatek_sygnal}")

    czas_poczatek_sygnal = time.time()
    widmo = widmo_amplitudowe(v)
    print(f"Czas DFT dla sygnału v: {time.time() - czas_poczatek_sygnal}")

    u = wykres_u(fs)
    czas_poczatek_sygnal = time.time()
    widmo = widmo_amplitudowe(u)
    print(f"Czas DFT dla sygnału u: {time.time() - czas_poczatek_sygnal}")

    wykres_b(wykresy=False, pomiar_czasu=True)

    czas_dft = time.time() - czas_poczatek

    print(f"\nSumaryczny czas FFT: {czas_fft}\n"
          f"Sumaryczny czas DFT: {czas_dft}")

zadanie_4()
