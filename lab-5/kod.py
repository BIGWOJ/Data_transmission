import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

#Data rozpoczęcia: 10.05
#Data końcowa: 24.05

def wykres(x, y, nazwa):
    plt.plot(x, y)
    plt.title(f"{nazwa}")
    plt.savefig(nazwa)
    plt.show()

def sygnal_odniesienia(f, t, A=1, fi=0):
    return A * np.sin(2*np.pi * f * t + fi)

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

def oblicz_h(sygnal_p):
    return max(sygnal_p) * 2/3

def c_na_bity(zmodulowany, ilosc_bitow, ask=False):

    wynik = []
    licznik_zero = 0
    licznik_jeden = 0
    ilosc_probek_bit = len(zmodulowany) // ilosc_bitow

    for i in range(ilosc_bitow):
        for j in range(ilosc_probek_bit*i, ilosc_probek_bit*(i+1)):
            if zmodulowany[j] == 0:
                licznik_zero += 1
            if zmodulowany[j] == 1:
                licznik_jeden += 1

        if not ask:
            if licznik_zero / 3 > licznik_jeden:
                wynik.append('0')
            else:
                wynik.append('1')

        else:
            if licznik_zero > licznik_jeden:
                wynik.append('0')
            else:
                wynik.append('1')

        licznik_zero = 0
        licznik_jeden = 0

    return wynik

def porownaj_bity(wejscie, odczytane, nazwa):

    print(f"====={nazwa}=====\n"
          f"Strumień bitowy wejściowy: {wejscie}\n"
          f"Strumień bitowy odczytany poprzez demodulację: {odczytane}")

    if np.array_equal(wejscie, odczytane):
        print("Odczytane bity są identyczne z wejściowym strumieniem bitowym.\n")

    plt.step(np.arange(len(wejscie)), wejscie, label=f"{nazwa} - wektor wejściowy", linewidth=5)
    plt.step(np.arange(len(odczytane)), odczytane, label=f"{nazwa} - odczytany wektor bitowy")
    plt.gca().invert_yaxis()
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15))
    plt.savefig(f"porownanie_{nazwa}")
    plt.show()

def demodulator_ask_psk(zmodulowany, sygnal_odniesienia, Tb, ilosc_bitow, psk=False):

    demodulowany = zmodulowany * sygnal_odniesienia

    if psk:
        wykres(np.arange(len(demodulowany)), demodulowany, "psk_z")

    else:
        wykres(np.arange(len(demodulowany)), demodulowany, "ask_z")

    ilosc_probek_bit = len(zmodulowany) // ilosc_bitow
    wyniki_calka = []
    ilosc_bitow = len(zmodulowany) // ilosc_probek_bit

    for i in range(ilosc_bitow):
        poczatek = i * ilosc_probek_bit
        koniec = poczatek + 2
        for _ in range(ilosc_probek_bit):
            calka = scipy.integrate.simpson(demodulowany[poczatek: koniec])
            koniec += 1
            wyniki_calka.append(calka)

    if psk:
        wykres(np.arange(len(wyniki_calka)), wyniki_calka, 'psk_p')

    else:
        wykres(np.arange(len(wyniki_calka)), wyniki_calka, 'ask_p')

    if psk:
        c_tab = np.where(np.array(wyniki_calka) < 0, 1, 0)
        wykres(np.arange(len(c_tab)), c_tab, "psk_c")

    else:
        h = oblicz_h(wyniki_calka)
        c_tab = np.where(np.array(wyniki_calka) >= h, 1, 0)
        wykres(np.arange(len(c_tab)), c_tab, "ask_c")

    odczytane = c_na_bity(c_tab, ilosc_bitow)

    return odczytane

def demodulator_fsk(zmodulowany, sygnal_odniesienia_1, sygnal_odniesienia_2, ilosc_bitow):

    x_1 = zmodulowany * sygnal_odniesienia_1
    x_2 = zmodulowany * sygnal_odniesienia_2
    wykres(np.arange(len(zmodulowany)), zmodulowany, "fsk_z")
    wykres(np.arange(len(x_1)), x_1, "fsk_x1(t)")
    wykres(np.arange(len(x_2)), x_2, "fsk_x2(t)")

    ilosc_probek_bit = len(zmodulowany) // ilosc_bitow
    wyniki_calka_x1 = []
    wyniki_calka_x2 = []
    ilosc_bitow = len(zmodulowany) // ilosc_probek_bit

    for i in range(ilosc_bitow):
        poczatek = i * ilosc_probek_bit
        koniec = poczatek + 2
        for _ in range(ilosc_probek_bit):
            calka_x1 = scipy.integrate.simpson(x_1[poczatek: koniec])
            calka_x2 = scipy.integrate.simpson(x_2[poczatek: koniec])
            koniec += 1
            wyniki_calka_x1.append(calka_x1)
            wyniki_calka_x2.append(calka_x2)

    wykres(np.arange(len(wyniki_calka_x1)), wyniki_calka_x1, 'fsk_p1')
    wykres(np.arange(len(wyniki_calka_x2)), wyniki_calka_x2, 'fsk_p2')

    p_tab = np.dot(-1, wyniki_calka_x1) + wyniki_calka_x2
    wykres(np.arange(len(p_tab)), p_tab, 'fsk_p')

    c_tab = np.where(p_tab > 0, 1, 0)
    wykres(np.arange(len(c_tab)), c_tab, "fsk_c")

    odczytane = c_na_bity(c_tab, ilosc_bitow, ask=False)

    return odczytane

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

    wzorcowy_sygnal = sygnal_odniesienia(fn, czestotliwosci_modulacja, A_1)

#ASK

    zmodulowany = ASK(A_1, A_2, strumien_bitowy, fn, czestotliwosci_modulacja, N)
    wykres(np.arange(len(wzorcowy_sygnal)), zmodulowany, 'ask_x')
    odczytane_bity = demodulator_ask_psk(zmodulowany, wzorcowy_sygnal, Tb, len(strumien_bitowy))
    porownaj_bity(strumien_bitowy, odczytane_bity, "ask")

#PSK

    zmodulowany = PSK(strumien_bitowy, fn, czestotliwosci_modulacja, N)
    wykres(np.arange(len(wzorcowy_sygnal)), zmodulowany, 'psk_x')
    odczytane_bity = demodulator_ask_psk(zmodulowany, wzorcowy_sygnal, Tb, len(strumien_bitowy), psk=True)
    porownaj_bity(strumien_bitowy, odczytane_bity, "psk")

#FSK

    fn_1 = (W + 1) / Tb
    fn_2 = (W + 2) / Tb

    zmodulowany = FSK(strumien_bitowy, fs, fn_1, fn_2, N)
    sygnal_odniesienia_1 = sygnal_odniesienia(fn_1, czestotliwosci_modulacja)
    sygnal_odniesienia_2 = sygnal_odniesienia(fn_2, czestotliwosci_modulacja)
    odczytane_bity = demodulator_fsk(zmodulowany, sygnal_odniesienia_1, sygnal_odniesienia_2, len(strumien_bitowy))
    porownaj_bity(strumien_bitowy, odczytane_bity, "fsk")

zadanie()
