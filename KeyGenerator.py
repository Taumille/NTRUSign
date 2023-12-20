import Polynome as pn


class KeyPair:
    def __init__(self, N, p, q):
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


if __name__ == "__main__":
    import time
    t = time.time()
    k = KeyPair(503, 2, (2, 8))
    print(f"Time to calculate key : {int((time.time()-t)*100)/100}s")
    print(f"Pub : ord({k.pub.ord()})")
    print(f"Priv : (ord({k.priv[0].ord()}), ord({k.priv[1].ord()}))")
