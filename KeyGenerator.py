import Polynome as pn


class KeyPair:
    def __init__(self, N, p, q, gen=False):
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


if __name__ == "__main__":
    import time
    t = time.time()
    k = KeyPair(503, 2, (2, 8), gen=True)
    print(f"Time to calculate key : {int((time.time()-t)*100)/100}s")
    k.export()
