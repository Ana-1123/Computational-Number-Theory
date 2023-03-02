from Cryptodome.Random import get_random_bytes
from Cryptodome.Util.number import getPrime, inverse
import random
import itertools as it
import numpy as np
import time


def init_param(m):
    # p = getPrime(161, randfunc=get_random_bytes)
    p = 11
    k = len(m) + 1
    return p, k


def nrortext_to_base(m, p):
    if m.isdigit():
        nr = int(m)
        res = list()
        while nr > 0:
            res.append(nr % p)
            nr //= p
        res = res[::-1]
    return res


def schemaHorner(m, n, x, p):
    m.insert(n + 1, 0)  # coeficientul a0=0
    value = m[0]  # ak-1
    for i in range(1, n + 1):
        value = value * x + m[i]
    return value % p


# y=[P(1),P(2),...,P(n)], n=k+2*s=k+2
def codificare(m, k, p, s=1):
    y = list()
    n = k + 2 * s
    for i in range(0, n):
        y.insert(i, schemaHorner(m, k - 1, i + 1, p))
    return y


def alterare(y, p):
    # modificam un element de pe o pozitie random
    err = random.randint(0, len(y) - 1)
    z = []
    for i in range(0, len(y)):
        if i == err:
            z.append((y[i] + 1) % p)
        else:
            z.append(y[i])
    return z, err


def coeficient_liberkk(z, k, p):
    start = time.time()
    Afrom = []
    for i in range(0, len(z)):
        Afrom.append(i + 1)  # Afrom {1,..n}
    a = []
    allposibleA = it.combinations(Afrom, k)
    for subset in list(allposibleA):
        subset = list(subset)
        fc = 0
        for i in subset:
            produs = 1
            for j in subset:
                if j != i:
                    produs *= (j * inverse(j - i, p))
            fc = (fc + z[i - 1] * produs % p) % p
        if fc == 0:
            a = subset.copy()
            break
    print("A", a)
    end = time.time()
    fc_time = end - start
    return a, fc_time


def coeficient_liberk(z, k, p):
    start = time.time()
    Afrom = []
    for i in range(0, len(z)):
        Afrom.append(i + 1)  # Afrom {1,..n}
    a = []
    allposibleA = it.combinations(Afrom, k)
    for subset in list(allposibleA):
        subset = list(subset)
        fc = 0
        for i in subset:
            produs = 1
            p_inv = 1
            for j in subset:
                if j != i:
                    produs *= j
                    p_inv *= (j - i)
            fc = (fc + z[i - 1] * (produs * inverse(p_inv, p)) % p) % p
        if fc == 0:
            a = subset.copy()
            break
    print("A", a)
    end = time.time()
    fc_time = end - start
    return a, fc_time


def coeficient_liber1(z, k, p):
    start = time.time()
    Afrom = []
    for i in range(0, len(z)):
        Afrom.append(i + 1)  # Afrom {1,..n}
    a = []
    fc = 0
    allposibleA = it.combinations(Afrom, k)
    for subset in list(allposibleA):
        subset = list(subset)
        fc = 0
        sum = 0
        for i in subset:
            produs_i = z[i - 1]
            p_j = 1
            for j in subset:
                if j != i:
                    produs_i *= (j - i)
                    p_j *= j
            sum += produs_i
        fc = p_j * inverse(sum, p)
        fc %= p
        print(fc)
        if fc == 0:
            a = subset.copy()
            break
    print("A", a)
    end = time.time()
    fc_time = end - start
    return a, fc_time


def comparare_calcul_fc(z, k, p):
    a, fc_kk = coeficient_liberkk(z, k, p)
    a, fc_k = coeficient_liberk(z, k, p)
    #a, fc_1 = coeficient_liber1(z, k, p)
    print("Calcul coeficientului liber varianta folosing k(k-1) inversari: ", fc_kk,
          "\nCalcul coeficientului liber varianta folosing k inversari: ", fc_k#,
          #"\nCalcul coeficientului liber varianta folosing 1 inversare: ", fc_1
          )


def multiply(x, y):
    x1 = np.array(x)
    x2 = np.array(y)
    mul = np.polymul(x2, x1)
    for i in range(0, len(mul)):
        mul[i] = mul[i] % p
    return mul


def add(A, B, m, n, p):
    size = max(m, n)
    sum = [0 for i in range(size)]
    for i in range(0, m, 1):
        sum[i] = A[i]
    for i in range(n):
        sum[i] += B[i]

    for i in range(0, len(sum)):
        sum[i] = sum[i] % p
    return sum


def decodificare(z, k, p):
    result = []
    a, fc_time = coeficient_liberk(z, k, p)
    for i in a:
        produs = [1]
        for j in a:
            if j != i:
                res = []
                inv = inverse(i - j, p)
                res.append(inv % p)  # (x-j)*inv(i-j)=x*inv(i-j)+(-j)*inv(i-j)
                res.append(((- j) * inv) % p)
                produs = multiply(produs, res)
        for l in range(0, len(produs)):
            produs[l] = produs[l] * z[i - 1] % p
        result = add(result, produs, len(result), len(produs), p)
    if (len(result) > 0):
        a0 = result.pop()
    print("m =", result)


if __name__ == '__main__':
    m = input("Introduceti un text:")
    p, k = init_param(m)
    vec = nrortext_to_base(m, p)
    print(vec)
    y = codificare(vec, k, p)
    print(codificare(vec, k, p))
    z, err = alterare(y, p)
    comparare_calcul_fc(z, k, p)
    decodificare(z, k, p)
