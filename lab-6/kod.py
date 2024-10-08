import random
import numpy as np

#Data rozpoczęcia: 20.05
#Końcowa data: 03.06

def koder_hamming_7_4(slowo):

    wynik = np.zeros(7)
    indeksy = [2,4,5,6]

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
        print(f"Zakodowane słowo przed korekcją: {zakodowane}")
        zakodowane[S-1] = (zakodowane[S-1] + 1) % 2
        print(f"Konieczność korekcji w miejscu indeksu: {S-1}\n"
              f"Zakodowane słowo po korekcji {zakodowane}")

    odkodowane = np.array([zakodowane[2], zakodowane[4], zakodowane[5], zakodowane[6]], dtype='int')

    return odkodowane

def koder_hamming_15_11(slowo):

    I = np.eye(11, dtype=int)
    #Macierz następnych liczb zapisanych binarnie bez wielokrotności 2
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

    return np.array(c), P

def dekoder_hamming_15_11(zakodowane, macierz_p):

    n = 15
    k = 11
    I = np.eye(n - k, dtype=int)
    H = np.hstack((I, macierz_p.transpose()))
    wektor_syndromu = np.dot(zakodowane, H.transpose()) % 2

    S = 0
    for i in range(len(wektor_syndromu)):
        S += wektor_syndromu[i] * 2 ** i

    if S != 0:
        print(f"Zakodowane słowo przed korekcją: {zakodowane}")
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
        print(f"Konieczność korekcji w miejscu indeksu: {korekta}\n"
              f"Zakodowane słowo po korekcji {zakodowane}")

    odkodowane = zakodowane[n-k:]

    return np.array(odkodowane)

def zadanie_1(slowo):

    print(f"\n==========Wersja 7_4==========\n"
          f"Słowo wejściowe: {slowo}\n")

    print(f"========== I przypadek - słowo bez zmian ==========")
    zakodowane = koder_hamming_7_4(slowo)
    odkodowane = dekoder_hamming_7_4(zakodowane)
    print(f"Słowo zakodowane: {zakodowane}\n"
          f"Słowo odkodowane: {odkodowane}")

    if np.array_equal(slowo, odkodowane):
        print("\n====WNIOSEK: Słowo odkodowane jest taki same jak oryginalne====\n")

    print(f"========== II przypadek - słowo z losową zmianą ==========")
    zanegowany_bit = random.randint(0,6)
    print(f"Losowo wybrany indeks, zanegowanego bitu w zakodowanym słowie: {zanegowany_bit}")
    zakodowane[zanegowany_bit] = (zakodowane[zanegowany_bit] + 1) % 2
    odkodowane = dekoder_hamming_7_4(zakodowane)
    print(f"Słowo odkodowane: {odkodowane}")

    if np.array_equal(slowo, odkodowane):
        print("\n====WNIOSEK: Słowo odkodowane jest taki same jak oryginalne====\n")

def zadanie_2(slowo):

    c, macierz_p = koder_hamming_15_11(slowo)
    odkodowane = dekoder_hamming_15_11(c, macierz_p)
    print(f"\n==========Wersja 15_11==========\n"
          f"Słowo wejściowe: {slowo}\n")
    print(f"========== I przypadek - słowo bez zmian ==========")
    print(f"Zakodowane słowo: {c}\n"
          f"Odkodowane słowo: {odkodowane}")

    if np.array_equal(slowo, odkodowane):
        print("\n====WNIOSEK: Słowo odkodowane jest taki same jak oryginalne====\n")

    print(f"========== II przypadek - słowo z losową zmianą ==========")
    zanegowany_bit = random.randint(0,14)
    print(f"Losowo wybrany indeks, zanegowanego bitu w zakodowanym słowie: {zanegowany_bit}")
    c[zanegowany_bit] = (c[zanegowany_bit] + 1) % 2
    odkodowane = dekoder_hamming_15_11(c, macierz_p)
    print(f"Słowo odkodowane: {odkodowane}")

    if np.array_equal(slowo, odkodowane):
        print("\n====WNIOSEK: Słowo odkodowane jest taki same jak oryginalne====")

slowo_4 = np.array(list('1101'), dtype=int)
slowo_11 = np.array(list('10101001010'), dtype=int)

zadanie_1(slowo_4)
zadanie_2(slowo_11)
