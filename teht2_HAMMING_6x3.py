# 6x3 Hamming virheenkorjaus
# Koodille annetaan 3-bittinen viestisana, johon generoidaan 1 virhe, jonka koodi tunnistaa ja korjaa


import numpy as np
import math
from random import randint


def identiteettimatriisi(n):
    identiteetti = np.zeros((n, n), dtype=int)
    for i in range(0, n):
        identiteetti[i][i] = 1
    return identiteetti


def tarkistusmatriisi(g):
    return np.hstack((np.transpose(G), identiteettimatriisi(3)))


def modulo(M):  # Modulo kaikille yli 1 arvoille
    if len(M.shape) == 1:  # 1 Ulotteinen matriisi
        for i in range(0, len(M)):
            if M[i] > 1:
                M[i] = M[i] % 2
    else:  # Yli 1 ulotteinen matriisi
        for j in range(0, M.shape[1]):
            for i in range(0, len(M[j])):
                if M[j][i] > 1:
                    M[j][i] = M[j][i] % 2
    return M


def virheGenerointi(r):
    virhemaara = randint(1, 2)
    virheet = np.zeros(6, dtype=int)
    for i in range(0, virhemaara):
        virhekohta = randint(0, len(r) - 1)
        virheet[virhekohta] = 1
    return modulo(r + virheet)


def hammingDistMin(kaikki, c):
    distance = [0] * len(kaikki)
    for r in range(0, len(kaikki)):
        for i in range(0, len(kaikki[r])):
            if kaikki[r][i] != c[i]:
                distance[r] = distance[r] + 1
    minDistance = len(c)
    for d in distance:
        if d < minDistance and d != 0:
            minDistance = d
    print(distance)
    return minDistance


# Lähetettävä viesti
m = np.array([1, 0, 0])

# Generointi matriisi
G = np.array([
    [1, 0, 1],
    [1, 1, 0],
    [0, 1, 1]])

# Tarkistusmatriisi
H = tarkistusmatriisi(G)

# Viimeistele generointimatriisi muotoon [I|G]
G = np.hstack((identiteettimatriisi(3), G))

# Lähetettävä koodi
c = m.dot(G)
modulo(c)
print("Koodisana c =", c)
# Vastaanotettu viesti
r = virheGenerointi(c)
print("Vastaanotettu r =", r)

HrT = modulo(H.dot(np.transpose(r)))

leveys = H.shape[1]
kandidaatit = np.zeros((leveys,), dtype=int)
syndromit = np.zeros((3,), dtype=int)
johtajat = np.zeros((3,), dtype=int)
s = 0  # kandidaattien generoimiseen käytettävä muuttuja

# Käydään kaikki kandidaatit läpi kunnes löydetään 2^n-m johtajaa
while johtajat.shape[0] < (2 ** (leveys - H.shape[0])):
    # Kandidaatti
    kandidaatti = np.zeros((leveys,), dtype=int)
    if math.floor(s / leveys) > 0:
        kandidaatti[math.floor(s / leveys) - 1] = 1
    if kandidaatti[s % leveys] == 1:
        s = s + 1
    kandidaatti[s % leveys] = 1
    s = s + 1
    # Tarkistetaan onko kandidaatti jo listassa
    for k in kandidaatit:
        if np.array_equal(k, kandidaatti):
            # Jos on listassa, niin skipataan se kokonaan
            continue

    # Lisätään kandidaatti kandidaattien listaan
    kandidaatit = np.vstack((kandidaatit, kandidaatti))
    # Syndromi
    syndromi = modulo(H.dot(kandidaatti))
    syndromit = np.vstack((syndromit, syndromi))
    # Johtaja
    onJoListassa = 0
    for j in johtajat:  # Tarkistetaan onko syndromia ollut aikaisemmin
        if np.array_equal(j, syndromi):
            onJoListassa = 1
    if onJoListassa == 0:  # Jos ei ole listassa niin lisätään johtajiin
        johtajat = np.vstack((johtajat, syndromi))

ctahti = np.zeros((6,), dtype=int)
if not np.array_equal(HrT, np.zeros((3,), dtype=int)):
    print("Viestissä on yksi tai useampi virhe, yritetään korjata.")
    for i in range(0, len(kandidaatit)):
        if np.array_equal(syndromit[i], HrT):
            ctahti = modulo(r + kandidaatit[i])
            break

print("Dekoodattu sana c* =", ctahti)
print("Oletettu sana m* =", ctahti[0:3])
