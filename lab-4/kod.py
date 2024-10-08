import numpy as np
import matplotlib.pyplot as plt
import math
#Data rozpoczęcia: 26.04
#Data końcowa: 10.05

def wykres(x, y, nazwa):
    plt.plot(x, y)
    plt.title(f"{nazwa}")
    plt.savefig(nazwa)
    plt.show()

def ASK(A_1, A_2, b, f, t, N):

    wynik = []
    aktualny_bit = 0
    liczba_probek_bit = N // len(b)
    aktualna_liczba_probek_bit = 0

    for i in range(N):
        if aktualna_liczba_probek_bit == liczba_probek_bit:
            aktualny_bit += 1
            aktualna_liczba_probek_bit = 0
        if b[aktualny_bit] == '0':
            wynik.append(A_1*np.sin(2*np.pi * f * t[i]))
        if b[aktualny_bit] == '1':
            wynik.append(A_2*np.sin(2*np.pi * f * t[i]))

        aktualna_liczba_probek_bit += 1

    return wynik

def PSK(b, f, t, N):

    wynik = []
    aktualny_bit = 0
    liczba_probek_bit = N // len(b)
    aktualna_liczba_probek_bit = 0

    for i in range(N):
        if aktualna_liczba_probek_bit == liczba_probek_bit:
            aktualny_bit += 1
            aktualna_liczba_probek_bit = 0
        if b[aktualny_bit] == '0':
            wynik.append(np.sin(2*np.pi * f * t[i]))
        if b[aktualny_bit] == '1':
            wynik.append(np.sin(2*np.pi * f * t[i] + np.pi))

        aktualna_liczba_probek_bit += 1

    return wynik

def FSK(b, fs, f_1, f_2, N):

    wynik = []
    aktualny_bit = 0
    liczba_probek_bit = N // len(b)
    aktualna_liczba_probek_bit = 0

    for i in range(N):
        if aktualna_liczba_probek_bit == liczba_probek_bit:
            aktualny_bit += 1
            aktualna_liczba_probek_bit = 0
        if b[aktualny_bit] == '0':
            wynik.append(np.sin(2 * np.pi * f_1 * i / fs))
        if b[aktualny_bit] == '1':
            wynik.append(np.sin(2 * np.pi * f_2 * i / fs))

        aktualna_liczba_probek_bit += 1

    return wynik

def szerokosc_pasma(decybele_widmo, decybele, czestotliwosci):

    poziom_odnesienia = max(decybele_widmo) - decybele
    obszar_szerokosci = np.where(decybele_widmo >= poziom_odnesienia)
    f_min = np.min(czestotliwosci[obszar_szerokosci])
    f_max = np.max(czestotliwosci[obszar_szerokosci])
    print(f"Szerokość pasma {(decybele)}dB: {f_max-f_min}")

def ascii_bity(tekst):

    wynik = []
    for litera in tekst:
        ascii = ord(litera) #Zapis litery do ASCII
        binarnie = bin(ascii)[2:] #[2:], żeby odciąc niepotrzebne 0b z zapisu binarneg
        if len(binarnie) < 7:
            roznica_dlugosci = 7-len(binarnie)
            for _ in range(roznica_dlugosci):
                wynik.append(0)
        wynik.extend(binarnie)

    return wynik

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

def zadanie():

    T = 1
    A_1 = 2
    A_2 = 5
    fs = 500
    #docelowa częstotliwość
    W = 2
    #ograniczona liczba bitów
    B = 10
    #czas pojedynczego bitu
    Tb = T / B
    fn = W * (1 / Tb)
    N = T * fs
    czestotliwosci_modulacja = np.arange(N) / fs
    strumien_bitowy = ascii_bity('przykladowy tekst')

    if len(strumien_bitowy) > B:
        strumien_bitowy = strumien_bitowy[:B]

    decybele = [3,6,12]

#ASK

    zmodulowany = ASK(A_1, A_2, strumien_bitowy, fn, czestotliwosci_modulacja, N)
    czas = np.linspace(0, T, len(zmodulowany))
    wykres(czas, zmodulowany, 'ASK')

    widmo = widmo_amplitudowe(zmodulowany)
    decybele_widmo = skala_decybelowa(widmo)
    czestotliwosci = np.arange(len(widmo))
    wykres(czestotliwosci, decybele_widmo, 'ASK - widmo')

    print("ASK")
    for decybel in decybele:
        szerokosc_pasma(decybele_widmo, decybel, czestotliwosci)

#PSK

    zmodulowany = PSK(strumien_bitowy, fn, czestotliwosci_modulacja, N)
    czas = np.linspace(0, T, len(zmodulowany))
    wykres(czas, zmodulowany, 'PSK')

    widmo = widmo_amplitudowe(zmodulowany)
    decybele_widmo = skala_decybelowa(widmo)
    czestotliwosci = np.arange(len(widmo))
    wykres(czestotliwosci, decybele_widmo, 'PSK - widmo')

    print('\nPSK')
    for decybel in decybele:
        szerokosc_pasma(decybele_widmo, decybel, czestotliwosci)


#FSK

    fn1 = (W + 1) / Tb
    fn2 = (W + 2) / Tb

    zmodulowany = FSK(strumien_bitowy, fs, fn1, fn2, N)
    czas = np.linspace(0, T, len(zmodulowany))
    wykres(czas, zmodulowany, "FSK")

    widmo = widmo_amplitudowe(zmodulowany)
    decybele_widmo = skala_decybelowa(widmo)
    czestotliwosci = np.arange(len(widmo))
    wykres(czestotliwosci, decybele_widmo, 'FSK - widmo')

    print('\nFSK')
    for decybel in decybele:
        szerokosc_pasma(decybele_widmo, decybel, czestotliwosci)


zadanie()
