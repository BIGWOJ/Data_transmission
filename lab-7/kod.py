import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
import binascii

def sygnal_odniesienia(f, t, A=1, fi=0):
    return A * np.sin(2*np.pi * f * t + fi)

def ascii_bity(tekst):
    wynik = []
    for litera in tekst:
        ascii = ord(litera)
        binarnie = bin(ascii)[2:].zfill(7)
        wynik.extend([int(bit) for bit in binarnie])
    return np.array(wynik, dtype=int)

def bity_ascii(ciag_binarny):
    byte_chunks = [ciag_binarny[i:i + 8] for i in range(0, len(ciag_binarny), 8)]
    ascii_string = ''.join([chr(int(''.join(map(str, byte)), 2)) for byte in byte_chunks])
    return ascii_string

def koder_hamming_7_4(slowo):
    wynik = np.zeros(7)
    indeksy = [2, 4, 5, 6]
    for i, indeks in enumerate(indeksy):
        wynik[indeks] = slowo[i]
    wynik[0] = (wynik[2] + wynik[4] + wynik[6]) % 2
    wynik[1] = (wynik[2] + wynik[5] + wynik[6]) % 2
    wynik[3] = (wynik[4] + wynik[5] + wynik[6]) % 2
    return np.array(wynik, dtype=int)

def dekoder_hamming_7_4(zakodowane):

    x1_prim = (zakodowane[2] + zakodowane[4] + zakodowane[6]) % 2
    x2_prim = (zakodowane[2] + zakodowane[5] + zakodowane[6]) % 2
    x4_prim = (zakodowane[4] + zakodowane[5] + zakodowane[6]) % 2
    x1 = (zakodowane[0] + x1_prim) % 2
    x2 = (zakodowane[1] + x2_prim) % 2
    x4 = (zakodowane[3] + x4_prim) % 2
    S = x1 * 2**0 + x2 * 2**1 + x4 * 2**2
    if S != 0:
        zakodowane[S-1] = (zakodowane[S-1] + 1) % 2
    odkodowane = np.array([zakodowane[2], zakodowane[4], zakodowane[5], zakodowane[6]], dtype='int')

    return odkodowane

def koder_hamming_15_11(slowo):

    I = np.eye(11, dtype=int)
    # Macierz następnych liczb zapisanych binarnie bez wielokrotności 2
    macierz_glowna = np.array([[0, 0, 1, 1],
                               [0, 1, 0, 1],
                               [0, 1, 1, 0],
                               [0, 1, 1, 1],
                               [1, 0, 0, 1],
                               [1, 0, 1, 0],
                               [1, 0, 1, 1],
                               [1, 1, 0, 0],
                               [1, 1, 0, 1],
                               [1, 1, 1, 0],
                               [1, 1, 1, 1]], dtype=int)

    P = np.array([macierz_glowna[:, -1], macierz_glowna[:, -2], macierz_glowna[:, 1], macierz_glowna[:, 0]]).transpose()
    G = np.hstack((P, I))
    c = np.dot(slowo, G) % 2

    return np.array(c)

def dekoder_hamming_15_11(zakodowane):

    n = 15
    k = 11
    I = np.eye(n - k, dtype=int)

    macierz_glowna = np.array([[0, 0, 1, 1],
                               [0, 1, 0, 1],
                               [0, 1, 1, 0],
                               [0, 1, 1, 1],
                               [1, 0, 0, 1],
                               [1, 0, 1, 0],
                               [1, 0, 1, 1],
                               [1, 1, 0, 0],
                               [1, 1, 0, 1],
                               [1, 1, 1, 0],
                               [1, 1, 1, 1]], dtype=int)

    macierz_p = np.array([macierz_glowna[:, -1], macierz_glowna[:, -2], macierz_glowna[:, 1], macierz_glowna[:, 0]]).transpose()

    H = np.hstack((I, macierz_p.transpose()))
    wektor_syndromu = np.dot(zakodowane, H.transpose()) % 2

    S = 0
    for i in range(len(wektor_syndromu)):
        S += wektor_syndromu[i] * 2 ** i

    if S != 0:
        indeksy = [3,5,6,7,8]
        if S in indeksy:
            if S == 3:
                zakodowane[4] = (zakodowane[4] + 1) % 2
                korekta = 4

            elif S == 8:
                zakodowane[3] = (zakodowane[3] + 1) % 2
                korekta = 3

            else:
                zakodowane[S] = (zakodowane[S] + 1) % 2
                korekta = S
        else:
            zakodowane[S-1] = (zakodowane[S-1] + 1) % 2
            korekta = S-1

    odkodowane = zakodowane[n-k:]

    return np.array(odkodowane)

def ASK(A_1, A_2, b, f, t, N):

    wynik = []
    liczba_probek_bit = N // len(b)
    for i in range(N):
        aktualny_bit = i // liczba_probek_bit
        if aktualny_bit >= len(b):
            break
        if b[aktualny_bit] == 0:
            wynik.append(A_1 * np.sin(2 * np.pi * f * t[i]))
        else:
            wynik.append(A_2 * np.sin(2 * np.pi * f * t[i]))

    return wynik

def PSK(b, f, t, N):

    print(len(b))
    wynik = []
    liczba_probek_bit = N // len(b)
    for i in range(N):
        aktualny_bit = i // liczba_probek_bit
        #print(aktualny_bit)
        if aktualny_bit >= len(b):
            break
        if b[aktualny_bit] == 0:
            wynik.append(np.sin(2 * np.pi * f * t[i]))
        else:
            wynik.append(np.sin(2 * np.pi * f * t[i] + np.pi))

    return wynik


    # wynik = []
    # liczba_probek_bit = N // len(b)
    # for i in range(len(b)):
    #     for _ in range(liczba_probek_bit):
    #         # if aktualny_bit >= len(b):
    #         #     break
    #         if b[i] == 0:
    #             wynik.append(np.sin(2 * np.pi * f * t[i]))
    #         else:
    #             wynik.append(np.sin(2 * np.pi * f * t[i] + np.pi))
    #
    # return wynik

def FSK(b, fs, f_1, f_2, N):

    wynik = []
    liczba_probek_bit = N // len(b)
    for i in range(N):
        aktualny_bit = i // liczba_probek_bit
        if aktualny_bit >= len(b):
            break
        if b[aktualny_bit] == 0:
            wynik.append(np.sin(2 * np.pi * f_1 * i / fs))
        else:
            wynik.append(np.sin(2 * np.pi * f_2 * i / fs))

    return wynik

def oblicz_h(sygnal_p):
    srednia = np.mean(sygnal_p)
    odchylenie_standardowe = np.std(sygnal_p)
    return srednia + odchylenie_standardowe

def bity_na_string(ciag_bitow):
    przedzial = [ciag_bitow[i:i + 8] for i in range(0, len(ciag_bitow), 8)]
    wynik = ''.join([chr(int(''.join(map(str, byte)), 2)) for byte in przedzial])
    return wynik

def normalizacja(sygnal):
    return sygnal / np.max(np.abs(sygnal))

def demodulator_ask_psk(zmodulowany, sygnal_odniesienia, ilosc_bitow, psk=False):

    zmodulowany = normalizacja(zmodulowany)
    demodulowany = zmodulowany * sygnal_odniesienia

    ilosc_probek_bit = len(zmodulowany) // ilosc_bitow
    wyniki_calka = []

    for i in range(ilosc_bitow):
        poczatek = i * ilosc_probek_bit
        koniec = poczatek + ilosc_probek_bit
        calka = scipy.integrate.simpson(demodulowany[poczatek: koniec])
        wyniki_calka.append(calka)

    if psk:
        c_tab = np.where(np.array(wyniki_calka) < 0, 1, 0)
    else:
        h = oblicz_h(wyniki_calka)
        c_tab = np.where(np.array(wyniki_calka) >= h, 1, 0)

    odczytane = np.array(c_tab, dtype=int)

    return odczytane

def demodulator_fsk(zmodulowany, sygnal_odniesienia_1, sygnal_odniesienia_2, ilosc_bitow):

    zmodulowany = normalizacja(zmodulowany)

    x_1 = zmodulowany * sygnal_odniesienia_1
    x_2 = zmodulowany * sygnal_odniesienia_2

    ilosc_probek_bit = len(zmodulowany) // ilosc_bitow
    wyniki_calka_x1 = []
    wyniki_calka_x2 = []

    for i in range(ilosc_bitow):
        poczatek = i * ilosc_probek_bit
        koniec = poczatek + ilosc_probek_bit
        calka_x1 = scipy.integrate.simpson(x_1[poczatek: koniec])
        calka_x2 = scipy.integrate.simpson(x_2[poczatek: koniec])
        wyniki_calka_x1.append(calka_x1)
        wyniki_calka_x2.append(calka_x2)

    p_tab = np.dot(-1, wyniki_calka_x1) + wyniki_calka_x2

    c_tab = np.where(p_tab > 0, 1, 0)

    odczytane = np.array(c_tab, dtype=int)

    return odczytane

def koder(hamming, slowo_informacyjne):

    #Podział na ramki
    dopelniona_ramka_len = 0
    if hamming == 1:
        ilosc_pelnych_ramek = len(slowo_informacyjne) // 4
        if ilosc_pelnych_ramek * 4 != len(slowo_informacyjne):
            dopelniona_ramka_len = len(slowo_informacyjne) - ilosc_pelnych_ramek * 4

    elif hamming == 2:
        ilosc_pelnych_ramek = len(slowo_informacyjne) // 11
        if ilosc_pelnych_ramek * 11 != len(slowo_informacyjne):
            dopelniona_ramka_len = len(slowo_informacyjne) - ilosc_pelnych_ramek * 11

    zakodowane = []

    #Kodowanie
    for i in range(ilosc_pelnych_ramek):
        if hamming == 1:
            zakodowane.append(koder_hamming_7_4(slowo_informacyjne[4*i : 4*(i+1)]))
        else:
            zakodowane.append(koder_hamming_15_11(slowo_informacyjne[11*i: 11*(i+1)]))

    if dopelniona_ramka_len != 0:
        dopelnienie = slowo_informacyjne[-dopelniona_ramka_len:]
        if hamming == 1:
            zera = np.array([0] * (4-dopelniona_ramka_len))
            dopelnienie = np.append(dopelnienie, zera)
            zakodowane.append(koder_hamming_7_4(dopelnienie))
        else:
            zera = np.array([0] * (11 - dopelniona_ramka_len))
            dopelnienie = np.append(dopelnienie, zera)
            zakodowane.append(koder_hamming_15_11(dopelnienie))

    zakodowane = np.concatenate(zakodowane)

    return zakodowane

def modulator(zakodowane, modulacja):

    print("Dobór wspólnych parametrów:")
    T = int(input("Czas T: "))
    fs = int(input("Częstotliwość fs: "))
    W = int(input("Parametr W: "))
    B = float(input("Czas pojedynczego bitu B: "))
    N = T * fs
    Tb = T / B
    fn = W * (1 / Tb)
    czestotliwosci_modulacja = np.arange(N) / fs

    fn1 = 0
    fn2 = 0
    if modulacja == 1:
        print("Dobór parametrów dla ASK:")
        A_1 = int(input("Częstotliwość A_1: "))
        A_2 = int(input("Częstotliwość A_2: "))

        zmodulowany = ASK(A_1, A_2, zakodowane, fn, czestotliwosci_modulacja, N)
        czas = np.linspace(0, T, len(zmodulowany))
        plt.plot(czas, zmodulowany)
        plt.title("Modulacja ASK")
        plt.show()

    if modulacja == 2:
        zmodulowany = PSK(zakodowane, fn, czestotliwosci_modulacja, N)
        czas = np.linspace(0, T, len(zmodulowany))
        plt.plot(czas, zmodulowany)
        plt.title("Modulacja PSK")
        plt.show()

    if modulacja == 3:
        fn1 = (W + 1) / Tb
        fn2 = (W + 2) / Tb
        zmodulowany = FSK(zakodowane, fs, fn1, fn2, N)
        czas = np.linspace(0, T, len(zmodulowany))
        plt.plot(czas, zmodulowany)
        plt.title("Modulacja FSK")
        plt.show()

    return zmodulowany, fs, fn, fn1, fn2, czestotliwosci_modulacja

def kanal_transmisyjny_szum(zmodulowany):

    srednia = 0
    odchylenie_standardowe = 1
    bialy_szum = normalizacja(np.random.normal(srednia, odchylenie_standardowe, size=len(zmodulowany)))
    alfa = float(input("Podaj wartość alfa [0,1], do wymnożenia przez biały szum: "))
    sygnal = zmodulowany + (bialy_szum * alfa)
    plt.plot(bialy_szum)
    plt.title("Biały szum")
    plt.show()
    plt.plot(sygnal)
    plt.title("Modulacja + biały szum")
    plt.show()

    return sygnal

def kanal_transmisyjny_tlumienie(zmodulowany, fs):

    beta = float(input("Podaj wartość beta [0,10], do wymnożenia w wykładniku podczas tłumienia sygnału: "))
    gt = np.exp(-beta*np.arange(len(zmodulowany)) / fs)
    sygnal = zmodulowany * gt

    plt.plot(gt)
    plt.title("gt")
    plt.show()

    plt.plot(sygnal)
    plt.title("Stłumiony sygnał")
    plt.show()

    return sygnal

def demodulator(modulacja, sygnal, fn, fn1, fn2, czestotliwosci_modulacja, zmodulowany, slowo_informacyjne):

    wzorcowy_sygnal = sygnal_odniesienia(fn, czestotliwosci_modulacja, max(zmodulowany))

    if modulacja == 1:
        odczytane_bity = demodulator_ask_psk(sygnal, wzorcowy_sygnal[:len(sygnal)], len(slowo_informacyjne))

    if modulacja == 2:
        odczytane_bity = demodulator_ask_psk(sygnal, wzorcowy_sygnal[:len(sygnal)], len(slowo_informacyjne), psk=True)

    if modulacja == 3:
        sygnal_odniesienia_1 = sygnal_odniesienia(fn1, czestotliwosci_modulacja)
        sygnal_odniesienia_2 = sygnal_odniesienia(fn2, czestotliwosci_modulacja)
        odczytane_bity = demodulator_fsk(sygnal,  sygnal_odniesienia_1[:len(sygnal)], sygnal_odniesienia_2[:len(sygnal)], len(slowo_informacyjne))

    return odczytane_bity

def dekoder(odczytane_bity, hamming):

    zdekodowane = []

    if hamming == 1:
        for i in range(len(odczytane_bity) // 7):
            zdekodowane.append(dekoder_hamming_7_4(odczytane_bity[7*i : 7*(i+1)]))

    else:
        for i in range(len(odczytane_bity) // 15):
            zdekodowane.append(dekoder_hamming_15_11(odczytane_bity[15*i : 15*(i+1)]))

    return np.concatenate(zdekodowane)

def zadanie(tlumienie=False, kaskadowo=False):

    tekst_wejsciowy = input("Wpisz teskt: ")
    slowo_informacyjne = ascii_bity(tekst_wejsciowy)
    plt.title("Binarny wektor wejściowy")
    plt.step(np.arange(len(slowo_informacyjne)), slowo_informacyjne)
    plt.show()

    hamming = int(input(("Wybór kodowania Hamminga:\n"
          "1. 7_4\n"
          "2. 15_11\n")))

    modulacja = int(input("Wybór modulacji:\n"
                          "1. ASK\n"
                          "2. PSK\n"
                          "3. FSK\n"))

    zakodowane = koder(hamming, slowo_informacyjne)
    zmodulowany, fs, fn, fn1, fn2, czestotliwosci_modulacja  = modulator(zakodowane, modulacja)

    if tlumienie:
        sygnal = kanal_transmisyjny_tlumienie(zmodulowany, fs)
        if kaskadowo:
            sygnal = kanal_transmisyjny_szum(sygnal)

    elif not tlumienie:
        sygnal = kanal_transmisyjny_szum(zmodulowany)
        if kaskadowo:
            sygnal = kanal_transmisyjny_tlumienie(sygnal, fs)

    zdemodulowany = demodulator(modulacja, sygnal, fn, fn1, fn2, czestotliwosci_modulacja, zmodulowany, slowo_informacyjne)
    odczytane_bity = dekoder(zdemodulowany, hamming)

    tekst_wyjsciowy = bity_na_string(odczytane_bity)

    blad = 0
    for i in range(len(odczytane_bity)):
        if odczytane_bity[i] != slowo_informacyjne[i]:
            blad+=1

    BER = blad / len(slowo_informacyjne)
    print(f"Tekst wyjściowy: {tekst_wyjsciowy}")
    print(f"BER: {BER}")

#zadanie()
#zadanie(tlumienie=True)
#zadanie(tlumienie=True, kaskadowo=True)
zadanie(kaskadowo=True)

#czas trwania bitu 1/4s