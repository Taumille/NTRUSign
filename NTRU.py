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


def Signing(k: KeyPair, D, N_bound):
    r = 0
    N = k.N
    q = k.q
    l_b = float('inf')
    while True:
        i = k.B
        m0 = H(D+r.to_bytes(10, 'big'), k.N)
        m = m0
        s = Polynome(N=N)
        si = Polynome(N=N)
        x = Polynome(N=N)
        y = Polynome(N=N)
        while i >= 1:
            # Perturb the point using the private lattice
            x.coeff = np.floor((m.star_multiply(k.priv[1][i])*(-1/q)).coeff)
            y.coeff = np.floor((m.star_multiply(k.priv[0][i])*(1/q)).coeff)

            si = x.star_multiply(k.priv[0][i]) + y.star_multiply(k.priv[1][i])
            m = si.star_multiply(k.priv[2][i] - k.priv[2][i-1]).mod(q)
            s = s + si
            i -= 1
        # Sign the perturbed point using the public lattice
        x.coeff = np.floor((m0.star_multiply(k.priv[1][0])*(-1/q)).coeff)
        y.coeff = np.floor((m0.star_multiply(k.priv[0][0])*(1/q)).coeff)
        s0 = x.star_multiply(k.priv[0][0]) + y.star_multiply(k.priv[1][0])
        s = s + s0

        # Check the signature
        b = NTRUNorm(s, s.star_multiply(k.pub) - m0, (0, q))
        if b < N_bound:
            break
        elif b < l_b:
            l_b = b
        r = r + 1

    return (D, r, s)


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
