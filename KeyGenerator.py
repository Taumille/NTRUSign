import Polynome as pn
import numpy as np


class KeyPair:
    def __init__(self, N=503, p=2, q=(2, 7), gen=False):
        if gen:
            while True:
                try:
                    f = pn.randomGenPoly(N=N)
                    bck = f.coeff
                    f.mod(p)
                    fp = f.inv(p)
                    print("Found fp")
                    f.coeff = bck
                    f.mod(q[0])
                    fq = f.inv(q[0], q[1])
                    print("Found fq")
                    break
                except Exception:
                    pass
            g = pn.randomGenPoly(N=N)

            h = g * fq * p
            h.mod(q[0]**q[1])

            self.pub = h
            self.priv = (f, fp)
        else:
            self.pub = None
            self.priv = None
        self.N = N
        self.p = p
        self.q = q

    def export(self, printk=True):
        if self.pub is None:
            print("No public key saved, please load or generate a public key")
            return
        s = "-----BEGIN NTRU PUBLIC KEY BLOCK-----\n\n"
        for c in self.pub.coeff:
            s += str(c) + "|"
        s = s[:-1]
        s += "\n=="+str(self.q)
        s += "\n\n-----END NTRU PUBLIC KEY BLOCK-----"

        if printk:
            print(s)
        return s

    def import_pub(self, s):
        cursor = 0
        for i in range(2):
            while (s[cursor] != '\n'):
                cursor += 1
            cursor += 1
        public_coeff = []
        sn = ""
        while s[cursor] != '\n':
            if s[cursor] == '|':
                public_coeff.append(int(sn))
                sn = ""
            else:
                sn += s[cursor]
            cursor += 1
        public_coeff.append(int(sn))
        sn = ""

        q = [None, None]
        cursor += 4
        while s[cursor] != ',':
            sn += s[cursor]
            cursor += 1
        q[0] = int(sn)
        sn = ""
        cursor += 1
        while s[cursor] != ')':
            sn += s[cursor]
            cursor += 1
        q[1] = int(sn)

        self.pub = pn.Polynome(N=len(public_coeff), q=q[0]**q[1])
        self.pub.coeff = np.array(public_coeff)
        self.q = tuple(q)
        self.N = self.pub.N


if __name__ == "__main__":
    import time
    t = time.time()
    k = KeyPair(503, 2, (2, 8), gen=True)
    print(f"Time to calculate key : {int((time.time()-t)*100)/100}s")
    s = k.export(False)
    k2 = KeyPair(gen=False)
    k2.import_pub(s)
    k2.export()
