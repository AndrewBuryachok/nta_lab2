def decomposition(n):
    print("decomposition")
    P = {}
    for p in range(2, int(n ** 0.5)):
        if n % p == 0:
            n //= p
            P[p] = 1
            while n % p == 0:
                n //= p
                P[p] += 1
    if n > 1:
        P[n] = 1
    print(P)
    return P

def build_table(al, n, P):
    print("table building")
    R = {}
    for p in P:
        R[p] = []
        for i in range(p):
            R[p].append(pow(al, (n - 1) * i // p, n))
    print(R)
    return R

def gcd(a, b):
    if a == 0:
        return (b, 0, 1)
    g, y, x = gcd(b % a, a)
    return (g, x - (b // a) * y, y)

def inverse(a, m):
    g, x, y = gcd(a, m)
    if g != 1:
        raise Exception("gcd != 1")
    a_1 = x % m
    print(str(a) + "^-1 mod " + str(m) + " = " + str(a_1))
    return a_1

def calc_xk(al_1, bt, n, p, Rp, X, k):
    for i in range(len(X)):
        bt = (bt * pow(al_1, X[i] * pow(p, i, n), n)) % n
    bt = pow(bt, n // (p ** (k + 1)), n)
    print("finding " + str(bt))
    for i in range(len(Rp)):
        if bt == Rp[i]:
            print("x" + str(k) + " = " + str(i))
            return i
    raise Exception("ERROR in calc xk")

def solve_system(system, mod):
    print("system solving")
    print(system)
    x = 0
    M = []
    i = 0
    for a, m in system:
        mi = mod // m
        mi_1 = inverse(mi, m)
        x = (x + a * mi * mi_1) % mod
        M.append((mi, mi_1))
        i += 1
    print(M)
    print("x = " + str(x))
    return x

def Silver_Pohlig_Hellman(al, bt, n):
    P = decomposition(n - 1)
    R = build_table(al, n, P)
    print("calculating xi")
    al_1 = inverse(al, n)
    system = []
    for p in P:
        X = []
        x = 0
        l = P[p]
        pl = p ** l
        print("p = " + str(p))
        for i in range(l):
            xi = calc_xk(al_1, bt, n, p, R[p], X, i)
            X.append(xi)
            x = (x + xi * p ** i) % pl
        system.append((x, pl))
        print("x = " + str(x) + " mod " + str(p))
    x = solve_system(system, n - 1)
    print(str(al) + "^" + str(x) + " = " + str(bt) + " mod " + str(n))
    print("result checking")
    print(pow(al, x, n) == bt)
    return x

#part 3

#Silver_Pohlig_Hellman(5, 11, 97)

#Silver_Pohlig_Hellman(3, 15, 43)

#Silver_Pohlig_Hellman(5, 11, 73)

#Silver_Pohlig_Hellman(28, 3, 43)

#Silver_Pohlig_Hellman(2, 1, 3)

#part 4

#Silver_Pohlig_Hellman(12, 13, 17)

#Silver_Pohlig_Hellman(139, 7, 157)

#Silver_Pohlig_Hellman(767, 33, 1109)

#Silver_Pohlig_Hellman(6785, 83, 10009)

#Silver_Pohlig_Hellman(65256, 17709, 100003)

#Silver_Pohlig_Hellman(30698, 29570, 104729)

#Silver_Pohlig_Hellman(893937, 902814, 1000999)

#Silver_Pohlig_Hellman(3462806, 7146766, 10005053)

#Silver_Pohlig_Hellman(3997814, 83667396, 100007377)

#Silver_Pohlig_Hellman(844235946, 316388906, 1000006099)

#Silver_Pohlig_Hellman(7328014316, 5693708660, 10000005553)

#Silver_Pohlig_Hellman(40493389860, 84077439232, 100000015277)
