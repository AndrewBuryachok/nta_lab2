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

def inverse(a, m):
    a_1 = pow(a, m - 2, m)
    print(str(a) + "^-1 mod " + str(m) + " = " + str(a_1))
    return a_1

def calc_xk(al, bt, n, p, Rp, X, k):
    left = bt
    al_1 = inverse(al, n)
    for i in range(len(X)):
        left = (left * pow(al_1, X[i] * pow(p, i, n), n)) % n
        #left = (left * pow(al, n - 1 - prime - X[i] * pow(p, i, n), n)) % n
    left = pow(left, n // (p ** (k + 1)), n)
    print("finding " + str(left))
    for i in range(len(Rp)):
        if left == Rp[i]:
            print("x" + str(k) + " = " + str(i))
            return i
    print("ERROR in calc xk")

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
    system = []
    for p in P:
        X = []
        x = 0
        l = P[p]
        pl = p ** l
        print("p = " + str(p))
        for i in range(l):
            xi = calc_xk(al, bt, n, p, R[p], X, i)
            X.append(xi)
            x = (x + xi * p ** i) % pl
        system.append((x, pl))
        print("x = " + str(x) + " mod " + str(p))
    x = solve_system(system, n - 1)
    print("result checking")
    print(pow(al, x, n) == bt)
    return x

Silver_Pohlig_Hellman(5, 11, 97)

Silver_Pohlig_Hellman(3, 15, 43)

Silver_Pohlig_Hellman(5, 11, 73)
