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


if __name__ == "__main__":
    infile = open("Alice.pdf", "rb")
    data = infile.read()
    convertToPolynome(data)
