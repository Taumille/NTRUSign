import random
import time
import numpy as np


class Polynome:
    def __init__(self, N=503):
        self.coeff = [0 for _ in range(0, N)]
        self.N = len(self.coeff)

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

    def __mul__(self, other):
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
    P = randomGenPoly(13, True, 8)
    P.inv(8)
