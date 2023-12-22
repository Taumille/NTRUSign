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

    def star_multiply(self, other):
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
        return self

    def evaluate(self, x):
        res = 0
        for i in range(len(self)):
            res += self.coeff[i] * x**i
        return res

    def inv(self, q, r=1):
        """
        Compute the inverse of self modulo the ideal (q^r, x^N-1)
        This algorithm is an application of the
        Extended Euclidean Algorithm
        """

        # Variable initialisation
        N = len(self)

        xp0 = Polynome(N=N)
        xp1 = Polynome(N=N)
        xp1.coeff[0] = 1
        yp0 = Polynome(N=N)
        yp0.coeff[0] = 1
        yp1 = Polynome(N=N)

        """
        We will compute xp*A + yp*B = 1
        As A is equivalent to 0 in this field
        We will have B*yp = 1 i.e. yp=B^-1
        """
        B = Polynome(N=N)
        B.coeff = self.coeff

        R = Polynome(N=N)
        R.coeff = self.coeff

        A = Polynome(N=N+1)
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
    while True:
        try:
            n = 16
            q = 5
            f = randomGenPoly(n, 7)
            g = randomGenPoly(n, 7)

            (F, G) = NTRUSolve(n, q, f, g)
            break
        except Exception:
            print("Exception")
    print(f'f={f}')
    print(f'g={g}')
    print(f'F={F}')
    print(f'G={G}')
    print(f'q = {modXnp1(f*G-g*F, n)} == {q}')
