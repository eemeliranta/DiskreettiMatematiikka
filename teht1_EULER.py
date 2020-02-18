# Eulerin funktio
# Laskee luvun suhteellisten alkulukujen määrän

from math import gcd


def eulerPhi(n):
    vastaus = 0
    for i in range(1, n):
        if (gcd(i, n) == 1):
            vastaus = vastaus + 1
    return vastaus


for i in range(2, 100):
    numero = i
    print("Luvun '" + str(numero) + "' suhteellisten alkulukujen määrä:", eulerPhi(numero))

n = 600

print(eulerPhi(n))
