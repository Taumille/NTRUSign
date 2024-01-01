from KeyGenerator import KeyPair
from Polynome import Polynome
import hashlib
import numpy as np


def H(s: bytes):
    """
    Convert the byte string to a polynomial
    using its sha1.
    """
    m = hashlib.sha1()
    m.update(s)
    m = m.hexdigest()
    N = 2**(int(np.ceil(np.log2(8*len(m)))))
    p = Polynome(N=N)
    m = ''.join(format(ord(x), 'b') for x in m)
    for i in range(len(m)):
        p.coeff[i] = int(m[i])
    return p


def NTRUNorm(P, Q, mod=(0, 0)):
    """
    Definition of a new norm which is the sum
    of the norm of the coefficient of the two
    polynomials
    """
    Pc, Qc = P.coeff, Q.coeff
    if mod[0] != 0 and mod[0]:
        q = abs(mod[0])
        Pc = Pc - q * (np.min(Pc) // q)
    if mod[1] != 0 and mod[1]:
        q = abs(mod[0])
        Qc = Qc - q * (np.min(Qc) // q)
    return np.norm(Pc) + np.norm(Qc)


if __name__ == "__main__":
    infile = open("Alice.pdf", "rb")
    data = infile.read()
    convertToPolynome(data)
