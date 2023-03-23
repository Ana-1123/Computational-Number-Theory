from Cryptodome.Random import get_random_bytes
from Cryptodome.Util.number import getPrime, inverse
import math
import time

list_lr_bin = []
list_lr_fixed_wind = []
list_lr_slid_wind = []


def generare_param():
    p = getPrime(700, randfunc=get_random_bytes)
    q = getPrime(700, randfunc=get_random_bytes)
    r = getPrime(700, randfunc=get_random_bytes)
    while p == q or q == r or p == r:
        p = getPrime(700, randfunc=get_random_bytes)
        q = getPrime(700, randfunc=get_random_bytes)
        r = getPrime(700, randfunc=get_random_bytes)
    # p = 3
    # q = 5
    # r = 7
    # n = 105
    # phi = 48
    # e = 11
    # d = 35
    # y = 17
    return p, q, r


def tcr_garner_3(p, q, r, x_p, x_q, x_r):
    c_2 = inverse(p, q)
    c_3 = (inverse(q, r) * inverse(p, r)) % r
    x = x_p
    u = (x_q - x) * c_2 % q
    x += u * p
    u = (x_r - x) * c_3 % r
    x += u * p * q
    return x


def tcr_garner_2(p2, q, x_p2, x_q):
    c_2 = inverse(p2, q)
    u = (x_q - x_p2) * c_2 % q
    x = x_p2 + u * p2
    return x


def multi_prime_rsa(p, q, r, e, d, y):
    start = time.time()
    x_p = pow(y % p, d % (p - 1), p)
    x_q = pow(y % q, d % (q - 1), q)
    x_r = pow(y % r, d % (r - 1), r)
    x = tcr_garner_3(p, q, r, x_p, x_q, x_r)
    print("Decriptare multiprime RSA(TCR cu Garner): ", x)
    end = time.time()
    fc_time = end - start
    print("Timpul decriptarii utilizand TCR cu algoritmul Garner:", fc_time)


def decriptare_librarii(y, n):
    start = time.time()
    decript = pow(y, d, n)
    print("Decriptare cu librarii: ", decript)
    end = time.time()
    fc_time = end - start
    print("Timpul decriptarii cu librarii:", fc_time)


def multi_power_rsa(p, q, r, e, d, y):
    start = time.time()
    x_q = pow(y % q, d % (q - 1), q)
    p2 = pow(p, 2)
    x0 = pow(y % p, d % (p - 1), p)
    epsilon = (y - (pow(x0, e, p2))) // p
    gamma = e * (pow(x0, e - 1, p2)) % p
    x1 = epsilon * inverse(gamma, p) % p
    x_p2 = x0 + p * x1
    x = tcr_garner_2(p2, q, x_p2, x_q)
    print("Decriptare multipower RSA(TCR cu Garner): ", x)
    end = time.time()
    fc_time = end - start
    print("Timpul decriptarii(multipower RSA) utilizand TCR cu algoritmul Garner:", fc_time)


def lr_bin_exp(x, n, m):
    n_b = list(map(int, bin(n).replace("0b", "")))
    k = len(n_b)
    y = 1
    lungimea = 0
    for i in range(0, k):
        y = (y * y) % m
        if n_b[i] == 1:
            y = (y * x) % m
            lungimea += 1
    print("Lungimea pentru left right binary: ", lungimea)
    list_lr_bin.append(lungimea)
    return y


def multi_prime_rsa_lr_bin(p, q, r, e, d, y):
    start = time.time()
    x_p = lr_bin_exp(y % p, d % (p - 1), p)
    x_q = lr_bin_exp(y % q, d % (q - 1), q)
    x_r = lr_bin_exp(y % r, d % (r - 1), r)
    x = tcr_garner_3(p, q, r, x_p, x_q, x_r)
    print("Decriptare multiprime RSA(TCR cu Garner) lr_binar: ", x)
    end = time.time()
    fc_time = end - start
    print("Timpul decriptarii utilizand TCR cu algoritmul Garner lr_binar:", fc_time)


def lr_fixed_wind(x, n, m, base):
    x_list = []
    for i in range(0, base):
        x_list.append(0)
    x_list[0] = 1
    lungime = 0
    for i in range(1, base):
        x_list[i] = (x_list[i - 1] * x) % m
        lungime += 1
    n = list(map(int, bin(n).replace("0b", "")))
    k = len(n)
    y = 1
    for i in range(k):
        y = pow(y, base, m)
        y = y * x_list[n[i]] % m
        lungime += 1
    print("Lungimea pentru left right fixed window: ", lungime)
    list_lr_fixed_wind.append(lungime)
    return y


def multi_prime_rsa_lr_fixed_wind(p, q, r, e, d, y):
    start = time.time()
    x_p = lr_fixed_wind(y % p, d % (p - 1), p, 2)
    x_q = lr_fixed_wind(y % q, d % (q - 1), q, 2)
    x_r = lr_fixed_wind(y % r, d % (r - 1), r, 2)
    x = tcr_garner_3(p, q, r, x_p, x_q, x_r)
    print("Decriptare multiprime RSA(TCR cu Garner) lr_fixed_wind: ", x)
    end = time.time()
    fc_time = end - start
    print("Timpul decriptarii utilizand TCR cu algoritmul Garner lr_fixed_wind:", fc_time)


def lr_slid_wind(x, n, m, w=4):
    xi = dict()
    xi[1] = x % m
    xi[2] = (xi.get(1) * xi.get(1)) % m
    lungimea = 1
    for omega in range(3, 2 ** w, 2):
        xi[omega] = xi.get(omega - 2) * xi.get(2) % m
        lungimea += 1
    n_b = list(map(int, bin(n).replace("0b", "")))
    n_b = list(reversed(n_b))
    y = 1
    i = len(n_b) - 1
    while i >= 0:
        if n_b[i] == 0:
            y = y * y % m
            i -= 1
        else:
            nr = 0
            j = max(i - w + 1, 0)
            while j < i and n_b[j] != 1:
                j += 1
            for lnr in range(1, i - j + 2):
                y = y * y % m
            for k in range(i, j - 1, -1):
                nr = nr * 10 + n_b[k]
            ni_nj = int(str(nr), 2)
            y = y * xi.get(ni_nj) % m
            lungimea += 1
            i = j - 1
    print("Lungimea pentru left-right sliding window: ", lungimea)
    list_lr_slid_wind.append(lungimea)
    return y


def multi_prime_rsa_lr_slid_wind(p, q, r, e, d, y):
    start = time.time()
    x_p = lr_slid_wind(y % p, d % (p - 1), p)
    x_q = lr_slid_wind(y % q, d % (q - 1), q)
    x_r = lr_slid_wind(y % r, d % (r - 1), r)
    x = tcr_garner_3(p, q, r, x_p, x_q, x_r)
    print("Decriptare multiprime RSA(TCR cu Garner) lr_slid_wind: ", x)
    end = time.time()
    fc_time = end - start
    print("Timpul decriptarii utilizand TCR cu algoritmul Garner lr_slid_wind:", fc_time)


if __name__ == '__main__':
    p, q, r = generare_param()
    e = getPrime(16, randfunc=get_random_bytes)
    # e = 11
    phi = (p - 1) * (q - 1) * (r - 1)
    while math.gcd(e, phi) != 1:
        e = getPrime(16, randfunc=get_random_bytes)
    n = p * q * r
    d = inverse(e, phi)
    x = input("Introduceti un numar:")
    # y = pow(int(x), e, n)
    # print("Numarul criptat: ", y)
    y = int(x)
    multi_prime_rsa(p, q, r, e, d, y)

    decriptare_librarii(y, n)

    multi_prime_rsa_lr_bin(p, q, r, e, d, y)

    multi_prime_rsa_lr_fixed_wind(p, q, r, e, d, y)

    multi_prime_rsa_lr_slid_wind(p, q, r, e, d, y)

    average_l_bin = sum(list_lr_bin) / len(list_lr_bin)
    average_l_fixed = sum(list_lr_fixed_wind) / len(list_lr_fixed_wind)
    average_l_slid = sum(list_lr_slid_wind) / len(list_lr_slid_wind)
    average_list = [average_l_bin, average_l_fixed, average_l_slid]
    denumiri_metode = ["Left-right binary method", "Left-right fixed window method", "Left-right sliding window method"]
    max_l = max(average_list)
    min_l = min(average_list)
    print("Lungime maxima are: ", denumiri_metode[average_list.index(max_l)])
    print("Lungime minima are: ", denumiri_metode[average_list.index(min_l)])
    n = p * p * q
    phi = (p * p - p) * (q - 1)
    while math.gcd(e, phi) != 1:
        e = getPrime(16, randfunc=get_random_bytes)
    d = inverse(e, phi)
    x = input("Introduceti un numar:")
    # y = pow(int(x), e, n)
    # print("Numarul criptat: ", y)
    y = int(x)
    multi_power_rsa(p, q, r, e, d, y)
    decriptare_librarii(y, n)
