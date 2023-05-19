from Cryptodome.Random import get_random_bytes
from Cryptodome.Util.number import getPrime, inverse
import random
import numpy as np
import math
from sympy import sieve


def Jacobi(a, n):
    t = 1
    while a != 0:
        while a % 2 == 0:
            a = a // 2
            if n % 8 == 3 or n % 8 == 5:
                t = -t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        a = a % n
    return t


def prime_factorization(nr):
    # nr = (nr1 ** e1) * (nr2 **e2) * ...
    divs = []
    for prime_div in sieve.primerange(2, nr):
        if nr % prime_div == 0:
            power = 0
            while nr % prime_div == 0:
                power += 1
                nr //= prime_div
            pair = (prime_div, power)
            divs.append(pair)
    return divs


def PrimitiveRootSafePrime1(p):
    a = random.randint(2, p - 1)
    if Jacobi(a, p) == -1:
        return a
    else:
        return p - a


def pow_mod(x, y, p):
    if y == 0:
        return 1
    res = 1
    x = x % p

    if (x == 0):
        return 0

    while (y > 0):
        if ((y & 1) == 1):
            res = (res * x) % p
        y = y >> 1
        x = (x * x) % p

    return res


def PrimitiveRoot1(factorization, p):
    # generam aleator alpha ∈ Z*(p):
    not_found = True
    while not_found:
        alpha = random.randrange(1, p)
        ok = 1
        for r in factorization:
            if pow_mod(alpha, (p - 1) // r[0], p) == 1:
                ok = 0
                break
        if ok == 1:
            not_found = False
    return alpha


def Shanks(a, b, p):
    m = int(np.ceil(math.sqrt(p - 1)))
    table = {}
    for j in range(0, m):
        table[j] = pow_mod(a, j, p)
    print("Lista de (j,a^j): ", table)
    table = dict(sorted(table.items(), key=lambda v: (v[1], v[0])))
    print("Lista de (j,a^j) ordonata: \n", table)
    am = pow_mod(inverse(a, p), m, p)  # a^-m
    y = b
    for i in range(0, m):
        for j in table.keys():
            if y == table[j]:
                return i * m + j
        y = y * am % p


def pollard(alpha, beta, p):
    xi = []
    ai = []
    bi = []
    x0 = 1
    xi.append(x0)
    a0 = 0
    ai.append(a0)
    b0 = 0
    bi.append(b0)
    xc1 = x0
    ac1 = a0
    bc1 = b0
    while True:
        if xc1 % 3 == 1:  # s1
            xc2 = (beta * xc1) % p
            ac2 = ac1 + 0
            bc2 = (bc1 + 1) % (p - 1)
        elif xc1 % 3 == 0:  # s2
            xc2 = (xc1 * xc1) % p
            ac2 = (2 * ac1) % (p - 1)
            bc2 = (2 * bc1) % (p - 1)
        else:  # s3
            xc2 = (alpha * xc1) % p
            ac2 = (ac1 + 1) % (p - 1)
            bc2 = bc1 + 0
        if xc2 in xi and xi.index(xc2) == len(xi)/2:
            print(xi)
            print(xi.index(xc2))
            i = xi.index(xc2)
            r = (bi[i] - bc2) % (p - 1)
            if r == 0:
                return "Failure"
            else:
                x = (inverse(r, (p - 1)) * (ac2 - ai[i])) % (p - 1)
                return x
        else:
            xi.append(xc2)
            ai.append(ac2)
            bi.append(bc2)
            xc1 = xc2 + 0
            ac1 = ac2 + 0
            bc1 = bc2 + 0


def TCRGaussAlgorithm(xi_list, pi_ei, n):
    result = 0
    for i in range(0, len(pi_ei)):
        ni = n // pi_ei[i]
        di = inverse(ni, pi_ei[i]) % pi_ei[i]
        result = result + xi_list[i] * ni * di % n
    return result


def SilverPohligHellman(factors, a, b, p):
    r = len(factors)
    xi_list = []
    for i in range(0, r):
        # (Compute xi = l0 + l1pi + ··· + lei−1pei−1 i, where xi = x mod pei
        q = factors[i][0]  # pi
        e = factors[i][1]  # ei
        l = []
        y = 1
        l.append(0)
        alpha = pow_mod(a, (p - 1) // q, p)
        for j in range(1, e + 1):
            if j == 1:
                y = y * pow_mod(a, l[j - 1] * inverse(q, p), p) % p
            else:
                y = y * pow_mod(a, l[j - 1] * pow(q, j - 2), p) % p
            beta = pow_mod(b * inverse(y, p), (p - 1) // pow(q, j), p)
            log = pollard(alpha, beta, p)
            l.append(log)
        suma = 0
        for index in range(1, e + 1):
            suma = suma + l[index] * pow(q, index - 1)
        xi_list.append(suma)
    pi_ei = []
    for i in range(0, len(factors)):  # n=p1^e1*p2^e2 ... x = xi mod pi^ei
        pi_ei.append(pow(factors[i][0], factors[i][1]))
    x = TCRGaussAlgorithm(xi_list, pi_ei, p - 1)
    return x


def countBits(number):
    return int((math.log(number) /
                math.log(2)) + 1)


def solovay_strassen(n, t=1000):
    for i in range(t):
        a = random.randint(2, n - 2)
        r = pow_mod(a, (n - 1) // 2, n)
        if r != 1 and r != n - 1:
            return False
        s = Jacobi(a, n)
        if r != s % n:
            return False
    return True


if __name__ == '__main__':
    # Shanks
    p = getPrime(16, randfunc=get_random_bytes)
    factorizationS = prime_factorization(p - 1)
    a = PrimitiveRoot1(factorizationS, p)
    b = random.randrange(1, p)
    print("Radacina primitiva", a)
    print("Skanks x=", Shanks(a, b, p))

    # Generare p
    primes_list = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    prim = 0
    while prim == 0:
        list_exponenti = []
        product = 1
        for i in range(0, 10):
            list_exponenti.append(random.randint(1, 10))
            product *= primes_list[i] ** list_exponenti[i]
        if solovay_strassen(product + 1) is True and countBits(product + 1) >= 128:
            prim = 1
    factorizationSPH = []
    for i in range(0, 10):
        pair = (primes_list[i], list_exponenti[i])
        factorizationSPH.append(pair)
    pSPH = product + 1

    print("Numar prim:", pSPH, "nr de biti: ", countBits(product + 1))
    bSPH = random.randrange(1, pSPH)
    aSPH = PrimitiveRoot1(factorizationSPH, pSPH)
    print("Radacina primitiva", aSPH)
    print("pollard", pollard(2, 11, 13))
    factorization = prime_factorization(250)
    print("SilverPohligHellman x=", SilverPohligHellman(factorization, 71, 210, 251))  # a,b,p
#  print("SilverPohligHellman x=", SilverPohligHellman(factorizationSPH, aSPH, bSPH, pSPH))
