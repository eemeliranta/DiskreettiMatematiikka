# BCH virheenkorjaus
# Koodille annetaan 7-bittinen viestisana, johon generoidaan 1 tai 2 virhettä, jotka koodi tunnistaa ja korjaa


import numpy as np
import math
from random import randint


# mod(2) polynomin kaikille arvoille
def modulo(p):
    m = [0] * len(p.c)
    for i in range(len(m)):
        m[i] = p.c[i] % 2
    return np.poly1d(m)


# Generoi satunnaisiin kohtiin 2 virhettä
def virheGenerointi(c):
    virheet = [0] * 15
    virhemaara = randint(1, 2)
    for i in range(0, virhemaara):
        virhekohta = randint(0, 14)
        virheet[virhekohta] = 1
    return modulo(np.polyadd(c, np.poly1d(virheet)))


# Purkaa yli x^3 arvot pienemmiksi, niitä vastaaviksi arvoiksi
def polyDeconstruct(p, g):
    out = np.poly1d([0])
    for i in range(p.o + 1):
        if p.c[p.o - i]:
            out = modulo(np.polyadd(out, g[i]))
    return out


# Kindof potenssilasku
def polyIncrease(p, p2):
    fullPoly = [0] * (15 - len(p.c))
    fullPoly2 = [0] * (15 - len(p2.c))
    positions = []
    positions2 = []
    for i in p.c:  # Full 15 length array
        fullPoly.append(int(i))
    for i in p2.c:  # Full 15 length array
        fullPoly2.append(int(i))

    for i in range(len(fullPoly)):  # Get the positions of ones
        if fullPoly[i] == 1:
            positions = positions + [14 - i]
    for i in range(len(fullPoly2)):  # Get the positions of ones
        if fullPoly2[i] == 1:
            positions2 = positions2 + [14 - i]

    for i in range(len(positions)):  # Add up the positions
        for j in range(len(positions2)):
            positions[i] = (positions[i] * positions2[j]) % 15

    fullPoly = [0] * 15
    for i in positions:  # Add ones to the positions
        fullPoly[14 - i] = 1
    return np.poly1d(fullPoly)


# Kindof jakolasku
def polyDecrease(p, p2):
    fullPoly = [0] * (15 - len(p.c))
    fullPoly2 = [0] * (15 - len(p2.c))
    positions = []
    positions2 = []
    for i in p.c:  # Full 15 length array
        fullPoly.append(int(i))
    for i in p2.c:  # Full 15 length array
        fullPoly2.append(int(i))

    for i in range(len(fullPoly)):  # Get the positions of ones
        if fullPoly[i] == 1:
            positions = positions + [14 - i]
    for i in range(len(fullPoly2)):  # Get the positions of ones
        if fullPoly2[i] == 1:
            positions2 = positions2 + [14 - i]

    for i in range(len(positions)):  # Add up the positions
        for j in range(len(positions2)):
            positions[i] = (positions[i] - positions2[j]) % 15

    fullPoly = [0] * 15
    for i in positions:  # Add ones to the positions
        fullPoly[14 - i] = 1
    return np.poly1d(fullPoly)


# Viestisana
m = np.poly1d([1, 1, 0, 0, 1, 0, 1])
print("Viestisana m \t\t\t\t=", m.c)

# Generointipolynomi
G = np.poly1d([1, 1, 1, 0, 1, 0, 0, 0, 1])

# Kaikki potenssien arvot
g = [np.poly1d([1]), np.poly1d([1, 0]), np.poly1d([1, 0, 0]), np.poly1d([1, 0, 0, 0]), np.poly1d([1, 1])]

# Generoi kaikkien potenssien arvot
for i in range(4, 14):  # Potenssien määrä
    a = np.poly1d(g[-1])  # aikaisempi polynomi
    a = np.polymul(a, np.poly1d([0, 0, 1, 0]))  # vastaa arvoa x
    if a.order == 4:  # Suurin polynomi sallittu
        a = modulo(np.polyadd(a, [-1, 0, 0, 1, 1]))  # Poistetaan suurin polynomi ja lisätään sitä vastaava arvo
    g.append(a)

# Koodisana
c = modulo(np.polymul(m, G))
# Vastaanotettu sana
r = virheGenerointi(c)
print("Koodisana c\t\t\t\t\t=", c.coefficients)
print("Vastaanotettu sana r\t\t=", r.c)

S1 = polyDeconstruct(r, g)
S3 = [0] * len(g)
for i in (range(len(r.c))):
    if (r.c[i]):
        # (x^3)^i kohtaan lisätään 1, vähän hankalasti koska lista on väärin päin
        pos = len(g) - 1 - ((r.o - i) * 3 % len(g))
        S3[pos] += 1
        S3[pos] = S3[pos] % 2  # Modulo 2 juuri lisätty arvo
S3 = polyDeconstruct(modulo(np.poly1d(S3)), g)

# Muutetaan polynomit 'puretusta' muodosta arvoon x^i (0 <= i <= 14)
for i in range(len(g)):
    if S1 == g[i]:
        S1 = [0] * 15
        S1[-i - 1] = 1
        S1 = np.poly1d(S1)
    if S3 == g[i]:
        S3 = [0] * 15
        S3[-i - 1] = 1
        S3 = np.poly1d(S3)

# S1 == S3 == 0
if S1 == np.poly1d([0]) and S3 == np.poly1d([0]):
    print("Vastaanotetussa sanassa ei ole virheitä.")
# S1^3 == S3
elif polyIncrease(S1, np.poly1d([1, 0, 0, 0])) == S3:
    print("Vastaanotetussa sanassa on yksi virhe.")
    # x^i = S1^i
    e = [0] * len(g)
    e[14 - S1.o] = 1
    e = np.poly1d(e)
else:
    print("Vastaanotetussa sanassa on kaksi virhettä.")
    if S3 == np.poly1d([0]):  # Erityis case jos S3 == 0
        nollaPisteet = [5, 10]
    else:
        beta = polyDecrease(S3, polyIncrease(S1, np.poly1d([1, 0, 0, 0])))  # beta = S3/S1^3
        # Kaikki mahdolliset nollapisteet
        npKaikki = [[1], [6, 13], [11, 12], [], [7, 9], [2, 8], [], [], [3, 14], [], [1, 4], [], [], [], []]
        nollaPisteet = npKaikki[beta.o]

    # S1 * x^i
    for i in range(len(nollaPisteet)):
        nollaPisteet[i] = (nollaPisteet[i] + S1.o) % len(g)  # nollapisteet * S1

    # e(x)
    e = [0] * len(g)
    for i in nollaPisteet:
        e[-i - 1] += 1
    e = np.poly1d(e)

ctahti = modulo(np.polyadd(r, e))
print("Korjattu koodisanassa c*\t=", ctahti.c)

mtahti = modulo(np.polydiv(ctahti, G)[0])

mtahtiOrig = mtahti
mtahti = []
for i in mtahtiOrig:
    mtahti.append(int(i))
mtahti = np.poly1d(mtahti)

print("Korjattu viestisana m* \t\t=", mtahti.c)
