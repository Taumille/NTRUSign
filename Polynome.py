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
            if self.coeff[i] != 0:
                s += str(self.coeff[i])+"X^"+str(i)+" + "
        s += str(self.coeff[0])
        return s

    def ord(self):
        for i in range(1, self.N+1):
            if self.coeff[self.N-i] != 0:
                return self.N-i
        return -1

    def mod(self, q):
        self.coeff %= q

    def inv(self, q):
        N = len(self)
        xp0 = Polynome(N=N, q=q)
        xp1 = Polynome(N=N, q=q)
        xp1.coeff[0] = 1
        yp0 = Polynome(N=N, q=q)
        yp0.coeff[0] = 1
        yp1 = Polynome(N=N, q=q)
        Rp = Polynome(N=N, q=q)
        Rp.coeff = self.coeff
        Pp = Polynome(N=len(P), q=q)
        Pp.coeff = self.coeff
        Qp = Polynome(N=N+1, q=q)
        Qp.coeff[N] = 1
        Qp.coeff[0] = -1 % q
        while np.linalg.norm(Rp.coeff) != 0:
            Dp, Rp = longDivide(Pp, Qp, q)
            (Pp, Qp) = (Qp, Rp)
            y0, y1 = yp1.coeff, (yp0 - Dp * yp1).coeff % q
            x0, x1 = xp1.coeff, (xp0 - Dp * xp1).coeff % q
            xp0.coeff = x0
            xp1.coeff = x1
            yp0.coeff = y0
            yp1.coeff = y1
        a = pow(int(Pp.coeff[0]), -1, q)
        if Pp.ord() > 0:
            raise Exception("Inversion Fails")
        res = yp0 * a
        res.N = N
        res.coeff = res.coeff[:N]
        res.coeff %= q
        return res


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

    Pinv = P.inv(3)
    print(Pinv)
    
    # print(f"Inverse mod 32 : {P.inv(32)}")
