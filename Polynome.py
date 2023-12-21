import numpy as np
import copy


class Polynome:
    def __init__(self, N=503, gen=False, o=0):
        self.coeff = np.array([0 for _ in range(0, N)])
        self.N = len(self.coeff)
        if gen:
            self.coeff[o] = 1

    def construct(self, coeff):
        self.coeff = np.array(coeff)
        self.N = len(self.coeff)

    def __len__(self):
        return self.N

    def __add__(self, other):
        res = Polynome(N=max(len(self), len(other)))

        for k in range(min(len(self), len(other))):
            res.coeff[k] = self.coeff[k] + other.coeff[k]
            res.coeff[k] = res.coeff[k]

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

    def __mul__(self, other):
        if isinstance(other, int) or isinstance(other, np.int64):
            res = Polynome(N=self.N)
            res.coeff = self.coeff * other
        elif isinstance(other, Polynome):
            res = Polynome(N=len(self) + len(other))
            for k in range(len(res)):
                for j in range(k+1):
                    try:
                        res.coeff[k] += self.coeff[j]*other.coeff[k-j]
                    except IndexError:
                        pass
        return res

    def star_multiply(self, other, q=-1):
        if q == -1:
            q = 2**32
        res = Polynome(N=max(len(self), len(other)))
        for k in range(len(res)):
            for i in range(len(res)):
                try:
                    if i <= k:
                        res.coeff[k] += self.coeff[i] * other.coeff[k-i]
                    else:
                        res.coeff[k] += self.coeff[i] * other.coeff[len(self)+k-i]
                    res.coeff[k] = res.coeff[k]
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

    def __truediv__(self, other):
        if isinstance(other, int):
            res = Polynome(N=self.N)
            res.coeff = self.coeff // other
            return res
        else:
            raise Exception(f"Can't divide polynome by {type(other)}")

    def ord(self):
        for i in range(1, self.N+1):
            if self.coeff[self.N-i] != 0:
                return self.N-i
        return -1

    def mod(self, q):
        self.coeff %= q

    def evaluate(self, x):
        res = 0
        for i in range(len(self)):
            res += self.coeff[i] * x**i
        return res


def xgcd(A, B):
    """
    Compute the extended euclidean algorithm
    on polynomial A and B of degree 0
    """
    if isinstance(A, Polynome) and isinstance(B, Polynome):
        a = A.coeff[0]
        b = B.coeff[0]
    else:
        a = A
        b = B

    if b == 0:
        return (1, 0, A)
    else:
        (x, y, g) = xgcd(b, a % b)
        return (y, x-(a//b)*y, g)


def modXnp1(f: Polynome, N):
    """
    Reduce f mod(X^n+1)
    """
    res = Polynome(N=N)
    for i in range(N):
        res.coeff[i] = f.coeff[i] - f.coeff[i+N]
    return res


def N(f: Polynome):
    """
    Compute the field norm of f
    """
    # f0 is the list of even coefficient
    f0 = Polynome(N=f.N//2)
    for i in range(f.N//2):
        f0.coeff[i] = f.coeff[2*i]
    # f1 is the list of even coefficient
    f1 = Polynome(N=f.N//2)
    for i in range(f.N//2):
        f1.coeff[i] = f.coeff[2*i+1]
    f02 = modXnp1(f0*f0, f0.N)
    f12 = modXnp1(f1*f1, f1.N)
    xf12 = f12 * Polynome(N=f.N, gen=True, o=1)
    Nf = f02 - xf12
    Nf.N = f.N//2
    Nf.coeff = Nf.coeff[:Nf.N]
    return Nf


def NTRUSolve(n, q, f, g):
    """
    Compute F and G satisfying f*F+g*G=q
    with f and g in Z[X]/(X^n+1)
    """
    if n == 1:
        (u, v, gcdfg) = xgcd(f, g)
        if gcdfg != 1:
            raise Exception(f"GCD(f,g) = {gcdfg} not equal 1")
            return
        Identity = Polynome(N=1, gen=True, o=0)
        (F, G) = (Identity*q*v, Identity*q*u)
        return (F, G)
    else:
        fp = N(f)
        gp = N(g)
        (Fp, Gp) = NTRUSolve(n//2, q, fp, gp)
        Fp2 = Polynome(N=Fp.N*2)
        for i in range(Fp.N):
            Fp2.coeff[i*2] = Fp.coeff[i]
            if i != 0:
                Fp2.coeff[i*2] *= -1
        F = modXnp1(Fp2*g, n)

        Gp2 = Polynome(N=Fp.N*2)
        for i in range(Gp.N):
            Gp2.coeff[i*2] = Gp.coeff[i]
        G = modXnp1(Gp2*f, n)
        #(F, G) = reduce(f, g, F, G)
        return (F, G)


def longDivide(A, B, q=503):
    # Compute the division A = QB+R
    Q = Polynome(N=len(A))
    R = copy.deepcopy(A)
    for i in range(A.ord()-B.ord(), -1, -1):
        try:
            Q.coeff[i] = R.coeff[B.ord()+i] * pow(int(B.coeff[B.ord()]), -1, q)
        except ValueError:
            raise Exception(f"Can't inverse {int(B.coeff[B.ord()])} in base {q}")
        for j in range(B.ord()+i, i-1, -1):
            R.coeff[j] = (R.coeff[j] - Q.coeff[i] * B.coeff[j-i]) % q

    return (Q, R)


def randomGenPoly(N=503, d=0):
    p = Polynome(N)
    ones_coeff = [k for k in range(N)]
    while len(ones_coeff) > d:
        ones_coeff.pop(np.random.randint(len(ones_coeff)))
    for i in ones_coeff:
        p.coeff[i] = 1
    return p


if __name__ == "__main__":
    P = Polynome(N=11)
    P.coeff = np.array([-1, 1, 1, 0, -1, 0, 1, 0, 0, 1, -1])
    P.mod(3)
    print(f'P : {P}')

    Pinv = P.inv(3)
    print(f'Inverse modulo 3 : {Pinv}')
    print(f'Test : p * Pinv = {P.star_multiply(Pinv, 3)}')
    P.coeff = np.array([-1, 1, 1, 0, -1, 0, 1, 0, 0, 1, -1])
    P.mod(2)
    print(f'P : {P}')
    Pinv = P.inv(2, 5)
    print(f'Inverse modulo 32 : {Pinv}')
    print(f'Test : p * Pinv = {P.star_multiply(Pinv, 32)}')
