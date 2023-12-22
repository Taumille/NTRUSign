import Polynome as pn
import numpy as np
from multiprocessing import Pool


def singleWorker(params):
    N, df, dg, q, t = params
    while True:
        try:
            ft = pn.randomGenPoly(N, df)
            gt = pn.randomGenPoly(N, dg)

            if t == 'transpose':
                f = ft
                fp = gt
            elif t == 'standard':
                (F, G) = pn.NTRUSolve(N, q, ft, gt)
                f = ft
                fp = F
            h = (f.inv(q).star_multiply(fp)).mod(q)
            break
        except Exception:
            pass
    return (f, fp, h)


class KeyPair:
    def __init__(self, N=256, df=128, dg=128, q=127, B=8, t='transpose', gen=False):
        if gen:
            f = [None for _ in range(B+1)]
            fp = [None for _ in range(B+1)]
            h = [None for _ in range(B+1)]

            params = [(N, df, dg, q, t) for _ in range(B+1)]
            with Pool(B) as p:
                res = p.map(singleWorker, params)

            f = [r[0] for r in res]
            fp = [r[1] for r in res]
            h = [r[2] for r in res]

            self.pub = h[0]
            self.priv = (f, fp, h)
        else:
            self.pub = None
            self.priv = None
        self.N = N
        self.B = B
        self.q = q
        self.df = df
        self.dg = dg

    def export(self, printk=True):
        if self.pub is None:
            print("No public key saved, please load or generate a public key")
            return
        s = "-----BEGIN NTRU PUBLIC KEY BLOCK-----\n\n"
        for c in self.pub.coeff:
            s += str(c) + "|"
        s = s[:-1]
        s += "\n=="+str(self.N)+"|"+str(self.B)+"|"+str(self.q)+"|"+str(self.df)+"|"+str(self.dg)
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
