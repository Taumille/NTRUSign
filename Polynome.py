import random
import numpy as np
import copy


class Polynome:
    def __init__(self, N=503, gen=False, o=0, q=2**32):
        self.coeff = np.array([0 for _ in range(0, N)])
        self.N = len(self.coeff)
        self.q = q
        if gen:
            self.coeff[o] = 1

    def construct(self, coeff):
        self.coeff = np.array(coeff)
        self.N = len(self.coeff)

    def __len__(self):
        return self.N

    def __add__(self, other):
        res = Polynome(N=max(len(self), len(other)), q=self.q)

        for k in range(min(len(self), len(other))):
            res.coeff[k] = self.coeff[k] + other.coeff[k]
            res.coeff[k] = res.coeff[k] % self.q

        if min(len(self), len(other)) == len(self):
            for k in range(len(other), len(self)):
                res.coeff[k] = other.coeff[k]
        else:
            for k in range(len(self), len(other)):
                res.coeff[k] = self.coeff[k]
        return res

    def __sub__(self, other):
        tmp = Polynome(N=other.N, q=self.q)
        tmp.coeff = -other.coeff % tmp.q
        return self + tmp

    def __mul__(self, other):
        if isinstance(other, int) or isinstance(other, np.int64):
            res = Polynome(N=self.N, q=self.q)
            res.coeff = self.coeff * other
            res.coeff % self.q
        elif isinstance(other, Polynome):
            res = Polynome(N=max(len(self), len(other)), q=self.q)
            for k in range(len(res)):
                for j in range(k+1):
                    try:
                        res.coeff[k] += self.coeff[j]*other.coeff[k-j]
                    except IndexError:
                        pass
        return res

    def star_multiply(self, other, q=-1):
        if q == -1:
            q = self.q
        res = Polynome(N=max(len(self), len(other)), q=q)
        for k in range(len(res)):
            for i in range(len(res)):
                try:
                    if i <= k:
                        res.coeff[k] += self.coeff[i] * other.coeff[k-i]
                    else:
                        res.coeff[k] += self.coeff[i] * other.coeff[len(self)+k-i]
                    res.coeff[k] = res.coeff[k] % q
                except IndexError:
                    res.coeff[k] += 0

        return res

    def __str__(self):
        s = ""
        for i in range(len(self)-1, 0, -1):
            s += str(self.coeff[i])+"X^"+str(i)+" + "
        s += str(self.coeff[0])
        return s

    def ord(self):
        for i in range(1, self.N+1):
            if self.coeff[self.N-i] != 0:
                return self.N-i
        return -1

    def mod(self, q):
        for k in range(self.ord() + 1):
            self.coeff[k] = self.coeff[k] % q

def addP(P,Q):
    n = deg(P)
    m = deg(Q)
    N = np.max([n,m])+1
    P2 = np.zeros(N)
    Q2 = np.zeros(N)
    P2[0:n+1] = P[0:n+1]
    Q2[0:m+1] = Q[0:m+1]
    return P2+Q2

def sousP(P,Q):
    return addP(P, (-1)*Q)

def prodP(P, Q): 
    n = deg(P)
    m = deg(Q)
    s = np.zeros(m+n+1)
    for i in range (n+1):
        for j in range (m+1):
            s[i+j] += P[i]*Q[j]
    return s

def estPasZero(P):
    return np.sum(P!=0)>0

def xgcd(a, b):
    """return (g, x, y) such that a*x + b*y = g = gcd(a, b)"""
    x0, x1, y0, y1 = 0, 1, 1, 0
    while a != 0:
        (q, a), b = divmod(b, a), a
        y0, y1 = y1, y0 - q * y1
        x0, x1 = x1, x0 - q * x1
    return b, x0, y0

def deg(P):
    if np.sum(P!=0):
        return np.max(np.nonzero(P))
    else:
        return 0

def divP (P, Q, p): 
    n = deg(P)
    m = deg(Q)
    D = np.zeros_like(P)
    Temp = np.zeros_like(P)
    while n>=m and estPasZero(P):
        t = P[n]*pow(int(Q[m]), -1, p) % p
        Temp[n-m] = t
        D = D + Temp # ce par quoi on divise
        P = sousP(P,prodP(Temp, Q)) %p
        Temp[n-m]=0
        n = deg(P)
    return D, P

def xgcdP(P, p):
    N = len(P)
    Q = np.zeros(N+1)
    Q[N] = 1
    Q[0] = -1
    Q = Q%p
    R = np.zeros_like(P)
    D = np.zeros_like(P)
    x0 = np.zeros_like(P)
    x1 = np.zeros_like(P)
    y0 = np.zeros_like(P)
    y1 = np.zeros_like(P)
    x1[0]=1
    y0[0]=1
    R = np.copy(P)
    while np.linalg.norm(R) != 0:
        D, R = divP(P, Q, p)
        P = np.copy(Q)
        Q = np.copy(R)
        y0, y1 = np.copy(y1), sousP(y0, prodP(D,y1)) %p
        x0, x1 = np.copy(x1), sousP(x0, prodP(D,x1)) % p
    a = pow(int(P[0]), -1, p)
    if deg(P) > 0:
        raise Exception("Inversion Fails")
    return a*y0 % p

def modinvP(P, N, p):
    Q = np.zeros(N+1)
    Q[N] = 1
    Q[0] = -1
    Q = Q%p
    return xgcdP(P, Q, p)

def longDivide(A, B, q=503):
    # Compute the division A = QB+R
    Q = Polynome(N=len(A), q=q)
    R = copy.deepcopy(A)
    for i in range(A.ord()-B.ord(), -1, -1):
        Q.coeff[i] = R.coeff[B.ord()+i] * pow(int(B.coeff[B.ord()]), -1, q)
        for j in range(B.ord()+i, i-1, -1):
            R.coeff[j] = (R.coeff[j] - Q.coeff[i] * B.coeff[j-i]) % q

    return (Q, R)


def randomGenPoly(N=503, inP=False, modq=2**32-1):
    p = Polynome(N, q=modq)
    for i in range(N):
        p.coeff[i] = random.randint(-1, 1)
    if inP:
        # P equiv Q[X^N-1]
        p.coeff[0] += p.coeff[-1]
        p.coeff[-1] = 0
    return p


if __name__ == "__main__":
    P = Polynome(N=11)
    P.coeff = np.array([-1, 1, 1, 0, -1, 0, 1, 0, 0, 1, -1])
    
    for k in range(len(P)):
        P.coeff %= 3
    print(P.coeff)

    Pinv = xgcdP(P.coeff, 3)
    print(Pinv.astype(int))
    
    # print(f"Inverse mod 32 : {P.inv(32)}")
