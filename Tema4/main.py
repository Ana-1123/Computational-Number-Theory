from Cryptodome.Random import get_random_bytes
from Cryptodome.Util.number import getPrime, inverse
import random
import numpy as np
import math
from sympy import sieve
from primePy import primes


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


def PrimitiveRoot1(factorization, p):
    # factorization = prime_factorization(p - 1)
    # generam aleator alpha ∈ Z*(p):
    not_found = True
    while not_found:
        alpha = random.randrange(1, p)
        ok = 1
        for r in factorization:
            if pow(alpha, (p - 1) // r[0], p) == 1:
                ok = 0
                break
        if ok == 1:
            not_found = False
    return alpha


def Skanks(a, b, p):
    m = int(np.ceil(math.sqrt(p - 1)))
    table = {}
    for j in range(0, m):
        table[j] = pow(a, j, p)
    print("Lista de (j,a^j): ", table)
    table = dict(sorted(table.items(), key=lambda v: (v[1], v[0])))
    print("Lista de (j,a^j) ordonata: \n", table)
    am = pow(inverse(a, p), m, p)  # a^-m
    y = b
    for i in range(0, m):
        for j in table.keys():
            if y == table[j]:
                return i * m + j
        y = y * am % p


def TCRGaussAlgorithm(xi_list, pi_ei, n):
    result = 0
    for i in range(0, len(pi_ei)):
        ni = n // pi_ei[i]
        di = inverse(ni, pi_ei[i]) % pi_ei[i]
        result = result + xi_list[i] * ni * di % n
    return result


def SilverPohligHellman(factors, a, b, p):
    # factors = prime_factorization(p - 1)
    r = len(factors)
    xi_list = []
    for i in range(0, r):
        # (Compute xi = l0 + l1pi + ··· + lei−1pei−1 i, where xi = x mod pei
        q = factors[i][0]  # pi
        e = factors[i][1]  # ei
        l = []
        y = 1
        l.append(0)
        alpha = pow(a, (p - 1) // q, p)
        for j in range(1, e + 1):
            if j == 1:
                y = y * pow(a, l[j - 1] * inverse(q, p), p) % p
            else:
                y = y * pow(a, l[j - 1] * pow(q, j - 2), p) % p
            beta = pow(b * inverse(y, p), (p - 1) // pow(q, j), p)
            log = Skanks(alpha, beta, p)
            l.append(log)
        suma = 0
        for index in range(1, e + 1):
            suma = suma + l[index] * pow(q, index - 1)
        xi_list.append(suma)
    pi_ei = []
    for i in range(0, len(factors)):  # n=p1^e1*p2^e2 ... x = xi mod pi^ei
        pi_ei.append(pow(factors[i][0], factors[i][1]))
    print("xi", xi_list)
    x = TCRGaussAlgorithm(xi_list, pi_ei, p - 1)
    return x


def countBits(number):
    return int((math.log(number) /
                math.log(2)) + 1)

def SolovayStrassen(n, t=100):
    # t>=1 ("security parameter") - nr iteratii
    for i in range(0, t):
        a = random.randint(2, n - 2)
        r = pow(a, int((n - 1) / 2), n)
        if r != 1 and r != n - 1:
            return 0
        s = Jacobi(a, n)
        if r != s % n:
            return 0
    return 1

# Test de primalitate optimizat
def is_prime(n):
    if n <= 3:
        return n > 1
    elif n % 2 == 0 or n % 3 == 0:
        return False

    for i in range(5, int(pow(n, 0.5)) + 1, 6):
        if n % i == 0 or n % (i + 2) == 0:
            return False

    return True

def get_large_prime():
    p = 16

    while not SolovayStrassen(p):
        p = 1
        for i in range(30):
            if is_prime(i):
                power = random.randint(1, 10)
                p *= i ** power
        p += 1

    return p


# 2,3,5,7,11,13,17,19
if __name__ == '__main__':
    # df = get_large_prime()
    # print(df)
    # print(countBits(df))

    p = getPrime(16, randfunc=get_random_bytes)
    factorizationS = prime_factorization(p - 1)
    a = PrimitiveRoot1(factorizationS, p)
    b = random.randrange(1, p)
    print("Radacina primitiva", a)
    print("Skanks x=", Skanks(a, b, p))
    # print("Skanks x=", Skanks(2, 11, 13))  # x=logab mod p
    # generam numere prime
    list_exponenti = []
    list_prime = [2, 3, 5, 7, 11, 13, 17, 19,23,29]
    for i in range(0, 10):
        list_exponenti.append(1)
    prim = 0
    sum = 1
    for i in range(0, 10):
        sum *= pow(list_prime[i], list_exponenti[i])
    while prim == 0:
        if primes.check(sum + 1) is True and countBits(sum) >= 64:
            prim = 1
        else:
            next_pow = random.randrange(0, 10)
            sum //= list_exponenti[next_pow] * list_prime[next_pow]
            list_exponenti[next_pow] = random.randint(1, 10)
            sum *= list_exponenti[next_pow] * list_prime[next_pow]
    pSPH = sum + 1
    print("Este prim:", pSPH)
    # factorizationSPH = []
    # for i in range(0, 10):
    #     pair = (list_prime[i], list_exponenti[i])
    #     factorizationSPH.append(pair)
    # primes_list = [2, 3, 5, 7, 11, 13, 17, 19]
    # prim = 0
    # while prim == 0:
    #     list_prime = []
    #     list_exponenti = []
    #     counter = random.randrange(1, 9)
    #     random_numbers = random.sample(range(0, 9), counter)
    #     for i in range(0, counter):
    #         list_prime.append(random_numbers[i])
    #         list_exponenti.append(1)
    #     sum = 1
    #     for i in range(0, counter):
    #         sum += pow(list_prime[i], list_exponenti[i])
    #     if primes.check(sum) is True and countBits(sum) >= 64:
    #         prim = 1
    #     else:
    #         next_pow = random.randrange(0, 29)
    #         sum -= list_exponenti[next_pow] * list_prime[next_pow]
    #         list_exponenti[next_pow] += 1
    #         sum += list_exponenti[next_pow] * list_prime[next_pow]
    # factorizationSPH = []
    # for i in range(0, 7):
    #     pair = (list_prime[i], list_exponenti[i])
    #     factorizationSPH.append(pair)
    # # print("Prim format: ", sum)
    # pSPH = sum
    # bSPH = random.randrange(1, sum)
    # aSPH = PrimitiveRoot1(factorizationSPH, sum)
    # # print("SilverPohligHellman x=", SilverPohligHellman(71, 210, 251))  # a,b,p
    # print("SilverPohligHellman x=", SilverPohligHellman(factorizationSPH, aSPH, bSPH, pSPH))
