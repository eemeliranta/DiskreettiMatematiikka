# 7x4 Hamming virheenkorjaus
# Koodille annetaan 4-bittinen viestisana, johon generoidaan 1 tai 2 virhettä,
# joita koodi yrittää korjata, mutta pystyy korjata oikein vain 1-bitin virheen


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


# Modulo kaikille yli 1 arvoille
def modulo(M):
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


# Generoi 1-2 virhettä annettuun viestiin
def virheGenerointi(r):
    virhemaara = randint(1, 2)
    virheet = np.zeros(7, dtype=int)
    for i in range(0, virhemaara):
        virhekohta = randint(0, len(r) - 1)
        virheet[virhekohta] = 1

    return modulo(r + virheet)


# Koodisanan minimietäisyys kaikista mahdollisista sanoista
def hammingDistMin(c):
    kaikki = list()
    # Generoi kaikki mahdolliset 4 bittiset sanat
    for i in range(0, 15):
        kaikki.append(list("".join(list(bin(i)[2:])).zfill(4)))
    # Muutetaan lista matriiseiksi
    kaikki = np.asarray(kaikki, dtype=np.int32)
    koodisanat = np.array([0, 0, 0, 0, 0, 0, 0])
    # Kaikki mahdolliset koodisanat
    for i in range(1, 15):
        koodisanat = np.vstack((koodisanat, modulo(kaikki[i].dot(G))))

    distances = [0] * len(koodisanat)
    # Käydään kaikki sanat läpi
    for r in range(0, len(koodisanat)):
        # Käydään kaikki sanan bitit läpi
        for i in range(0, len(koodisanat[r])):
            if koodisanat[r][i] != c[i]:
                distances[r] = distances[r] + 1
    minDistance = len(c)
    for d in distances:
        if d < minDistance and d != 0:
            minDistance = d
    print(distances)
    return minDistance


# Kysy viesti
m = np.array([1, 1, 0, 0])

# Generointi matriisi
G = np.array([
    [0, 1, 1],
    [1, 0, 1],
    [1, 1, 0],
    [1, 1, 1]])

# Tarkistusmatriisi
H = tarkistusmatriisi(G)

# Viimeistele generointimatriisi muotoon [I|G]
G = np.hstack((identiteettimatriisi(4), G))

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
x = 0

# Käydään kaikki kandidaatit läpi kunnes löydetään 2^n-m johtajaa
while johtajat.shape[0] < (2 ** (leveys - H.shape[0])):
    x = x + 1
    if x > 60:
        break

    # Kandidaatti
    kandidaatti = np.zeros((leveys,), dtype=int)

    if math.floor(s / leveys) > leveys:
        kandidaatti[math.floor(s / leveys ** 2) - 1] = 1  # 3 painoinen kandidaatti
    if kandidaatti[(math.floor(s / leveys) - 1) % leveys] == 1:
        s = s + 1

    if math.floor(s / leveys) > 0:
        kandidaatti[(math.floor(s / leveys) - 1) % leveys] = 1  # 2 painoinen kandidaatti
    if kandidaatti[s % leveys] == 1:
        s = s + 1

    kandidaatti[s % leveys] = 1  # 1 painoinen kandidaatti
    s = s + 1

    for k in kandidaatit:  # Tarkistetaan onko kandidaatti jo listassa
        if np.array_equal(k, kandidaatti):
            continue  # Jos on listassa, niin skipataan se kokonaan

    kandidaatit = np.vstack((kandidaatit, kandidaatti))  # Lisätään kandidaatti kandidaattien listaan

    # Syndromi
    syndromi = modulo(H.dot(kandidaatti))
    syndromit = np.vstack((syndromit, syndromi))

    # Johtaja
    onJoListassa = 0
    # Tarkistetaan onko syndromia ollut aikaisemmin
    for j in johtajat:
        if np.array_equal(j, syndromi):
            onJoListassa = 1

    # Jos ei ole listassa niin lisätään johtajiin
    if onJoListassa == 0:
        johtajat = np.vstack((johtajat, syndromi))

ctahti = np.zeros((6,), dtype=int)
if not np.array_equal(HrT, np.zeros((3,), dtype=int)):
    print("Viestissä on yksi tai useampi virhe, yritetään korjata.")
    for i in range(0, len(kandidaatit)):
        if np.array_equal(syndromit[i], HrT):
            ctahti = modulo(r + kandidaatit[i])
            break
print("Dekoodattu sana c* =", ctahti)
print("Oletettu sana m* =", ctahti[0:4])
# hammingDistMin(r)
