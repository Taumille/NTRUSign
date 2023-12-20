import random
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
        res = Polynome(N=max(len(self), len(other)), q=self.q)

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
        tmp = Polynome(N=other.N, q=self.q)
        tmp.coeff = -other.coeff
        return self + tmp

    def __mul__(self, other):
        if isinstance(other, int) or isinstance(other, np.int64):
            res = Polynome(N=self.N, q=self.q)
            res.coeff = self.coeff * other
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

    def inv(self, q, r=1):
        """
        Compute the inverse of self modulo the ideal (q^r, x^N-1)
        This algorithm is an application of the
        Extended Euclidean Algorithm
        """

        # Variable initialisation
        N = len(self)

        xp0 = Polynome(N=N, q=q)
        xp1 = Polynome(N=N, q=q)
        xp1.coeff[0] = 1
        yp0 = Polynome(N=N, q=q)
        yp0.coeff[0] = 1
        yp1 = Polynome(N=N, q=q)

        """
        We will compute xp*A + yp*B = 1
        As A is equivalent to 0 in this field
        We will have B*yp = 1 i.e. yp=B^-1
        """
        B = Polynome(N=N, q=q)
        B.coeff = self.coeff

        R = Polynome(N=N, q=q)
        R.coeff = self.coeff

        A = Polynome(N=N+1, q=q)
        A.coeff[N] = 1
        A.coeff[0] = -1 % q

        while np.linalg.norm(R.coeff) != 0:
            # Process to polynomial Long Division B/A
            Q, R = longDivide(B, A, q)
            (B, A) = (A, R)
            # Increment xp and yp according to EED
            yp0.coeff, yp1.coeff = yp1.coeff, (yp0 - Q * yp1).coeff % q
            xp0.coeff, xp1.coeff = xp1.coeff, (xp0 - Q * xp1).coeff % q

        if B.ord() > 0:
            raise Exception("Inversion Fails")
        # Format the result
        c = yp0 * pow(int(B.coeff[0]), -1, q)
        c.mod(q)
        # Truncate the result to ensure a length of N
        c.N = N
        c.coeff = c.coeff[:N]

        """
        If r>1 i.e. we're trying to calculate the inverse modulo p^r
        with p prime and r an integer.
        We use the strategy decribed on page 27 of "NTRU: A NEW HIGH SPEED
        PUBLIC KEY CRYPTOSYSTEM" by Jeffrey HOPSTEIN, Jill PIPHER and
        Joseph H. SILVERMAN.
        """
        if r != 1:
            p = q
            Identity = Polynome(N=N, gen=True, o=0)
            while p < q**r:
                self.q = p**2
                cp = (self.star_multiply(c, p**2) - Identity) / p
                c.q = p**2
                c = c - c.star_multiply(cp, p) * p
                p = p**2
            c.mod(q**r)
        return c


def longDivide(A, B, q=503):
    # Compute the division A = QB+R
    Q = Polynome(N=len(A), q=q)
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
    p = Polynome(N, q=2)
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
