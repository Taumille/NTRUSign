import random
import time
import numpy as np


class Polynome:
    def __init__(self, N=503, gen=False, o=0):
        self.coeff = np.array([0 for _ in range(0, N)])
        self.N = len(self.coeff)
        if gen:
            self.coeff[o] = 1

    def __len__(self):
        return self.N

    def __add__(self, other):
        res = Polynome(N=max(len(self), len(other)))

        for k in range(min(len(self), len(other))):
            res.coeff[k] = self.coeff[k] + other.coeff[k]

        if min(len(self), len(other)) == len(self):
            for k in range(len(other), len(self)):
                res.coeff[k] = other.coeff[k]
        else:
            for k in range(len(self), len(other)):
                res.coeff[k] = self.coeff[k]
        return res

    def __sub__(self, other):
        tmp = Polynome(N=other.N)
        tmp.coeff = -other.coeff
        return self + tmp

    def conv(self, other):
        res = Polynome(N=max(len(self), len(other)))
        for k in range(len(res)):
            for i in range(len(res)):
                if i <= k:
                    res.coeff[k] += self.coeff[i] * other.coeff[k-i]
                else:
                    res.coeff[k] += self.coeff[i] * other.coeff[len(self)+k-i]

        return res

    def __str__(self):
        s = ""
        for i in range(len(self)-1, 0, -1):
            s += str(self.coeff[i])+"X^"+str(i)+" + "
        s += str(self.coeff[0])
        return s

    def __mul__(self, other):
        if isinstance(other, int):
            res = Polynome(N=self.N)
            res.coeff = self.coeff * other

        elif isinstance(other, Polynome):
            res = Polynome(N=max(len(self), len(other)))
            for k in range(len(res)):
                for i in range(len(res)):
                    if i <= k:
                        res.coeff[k] += self.coeff[i] * other.coeff[k-i]
        else:
            raise TypeError(f"Unsupported operand type for Polynome : {type(other)}")

        return res


def euclidDiv(A_in, B_in, q=2**32):
    # Compute A//B with A and B two polynomes
    (A, B) = (A_in, B_in)
    Q = Polynome(A.N)
    R = A

    order = A.ord() - B.ord()

    for o in range(order, -1, -1):
        m = R.coeff[B.ord()+o]/B.coeff[B.ord()]
        if m < 0:
            m = np.ceil(m)
        else:
            m = np.floor(m)
        m = int(m) % q
        print(f"m : {m}")
        R = R - Polynome(N=A.N, gen=True, o=o) * B * m
        Q.coeff[o] = m

    return (R, Q)


def randomGenPoly(N=503, inP=False, modq=2**32-1):
    p = Polynome(N)
    for i in range(N):
        p.coeff[i] = random.randint(0, N)
    if inP:
        # P equiv Q[X^N-1]
        p.coeff[0] += p.coeff[-1]
        p.coeff[-1] = 0
    return p


if __name__ == "__main__":
    P = Polynome(7)
    P.coeff[5] = 2
    P.coeff[4] = 3
    P.coeff[1] = 2
    P.coeff[0] = 4
    S = Polynome(7)
    S.coeff[0] = 1
    S.coeff[2] = 3

    (R, Q) = euclidDiv(P, S)
    print(f"Q : {Q}")
    print(f"R : {R}")
