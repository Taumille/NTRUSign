from KeyGenerator import KeyPair
from Polynome import Polynome
import hashlib
import numpy as np


def H(s: bytes, N: int):
    """
    Convert the byte string to a polynomial
    using its sha1.
    """
    h = hashlib.sha1()
    i = 0
    m = ""
    while len(m) < N:
        h.update(s+str(i).encode("ascii"))
        m += h.hexdigest()
    p = Polynome(N=N)
    for i in range(len(m)):
        p.coeff[i % N] += ord(m[i])
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
        Pc = Pc % q
    if mod[1] != 0 and mod[1]:
        q = abs(mod[1])
        Qc = Qc % q
    res_p = np.sqrt(np.sum(np.square(Pc)) - np.square(np.sum(Pc))/len(Pc))
    res_q = np.sqrt(np.sum(np.square(Qc)) - np.square(np.sum(Qc))/len(Qc))
    return np.sqrt(res_p**2 + res_q**2)


if __name__ == "__main__":
    infile = open("Alice.pdf", "rb")
    data = infile.read()

    # Load key
    k = KeyPair()
    infile = open("key_pub.asc", "r")
    s = infile.read()
    infile = open("key_priv.asc", "r")
    s_priv = infile.read()
    k.import_pub(s)
    k.import_priv(s_priv)

    print("Keys imported")

    print(k.priv)
    print(Signing(k=k, D=data, N_bound=310))
