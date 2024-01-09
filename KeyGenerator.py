import Polynome as pn
import numpy as np
from multiprocessing import Pool


def singleWorker(params):
    """
    A single function that can be executed in parallel to accelerate
    Key creation
    """
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
    def __init__(self,
                 N=256,
                 q=127,
                 df=128,
                 dg=128,
                 B=8,
                 t='transpose',
                 gen=False,
                 name="User Name",
                 email="user@example.com"):
        """
        Create a key with the parameter passed to the constructor
        """
        if gen:
            f = [None for _ in range(B+1)]
            fp = [None for _ in range(B+1)]
            h = [None for _ in range(B+1)]

            nproc = max(24, B+1)

            params = [(N, df, dg, q, t) for _ in range(nproc)]
            i = 0
            res = []
            with Pool(nproc) as p:
                print(f"Launching {nproc} processes")
                for i in p.imap_unordered(singleWorker, params, chunksize=1):
                    res += [i]
                    if len(res) >= B:
                        break

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
        self.name = name
        self.email = email

    def export_pub(self, printk=True):
        """
        Export the key in a readable format either by returning it to the
        standard output or saving it in a string.
        """
        if self.pub is None:
            print("No public key saved, please load or generate a public key")
            return
        s = "-----BEGIN NTRU PUBLIC KEY BLOCK-----\n"
        s += self.name+"<"+self.email+">\n\n"
        for c in self.pub.coeff:
            s += str(c) + "|"
        s = s[:-1]
        s += "\n=="+str(self.N)+"|"+str(self.B)+"|"+str(self.q)+"|"+str(self.df)+"|"+str(self.dg)
        s += "\n\n-----END NTRU PUBLIC KEY BLOCK-----"

        if printk:
            print(s)
        return s

    def import_pub(self, s):
        """
        Import a key previously exported by the export method.
        """
        self.name = ""
        self.email = ""
        cursor = 0
        while s[cursor] != '\n':
            cursor += 1
        cursor += 1
        while s[cursor] != '<':
            self.name += s[cursor]
            cursor += 1
        cursor += 1
        while s[cursor] != '>':
            self.email += s[cursor]
            cursor += 1
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

        cursor += 3
        while s[cursor] != '|':
            sn += s[cursor]
            cursor += 1
        self.N = int(sn)
        sn = ""
        cursor += 1
        while s[cursor] != '|':
            sn += s[cursor]
            cursor += 1
        self.B = int(sn)
        sn = ""
        cursor += 1
        while s[cursor] != '|':
            sn += s[cursor]
            cursor += 1
        self.q = int(sn)
        sn = ""
        cursor += 1
        while s[cursor] != '|':
            sn += s[cursor]
            cursor += 1
        self.df = int(sn)
        sn = ""
        cursor += 1
        while s[cursor] != '\n':
            sn += s[cursor]
            cursor += 1
        self.dg = int(sn)

        self.pub = pn.Polynome(N=self.N)
        self.pub.coeff = np.array(public_coeff)

    def export_priv(self, printk=True):
        if self.priv is None:
            print("No priv key saved, please load or generate a public key")
            return
        s = "-----BEGIN NTRU PRIVATE KEY BLOCK-----\n"
        s += self.name+"<"+self.email+">\n\n"
        for i in range(self.B):
            for c in self.priv[0][i].coeff:
                s += str(c)+"|"
            s = s[:-1]
            s += '\n'

            for c in self.priv[1][i].coeff:
                s += str(c)+"|"
            s = s[:-1]
            s += '\n'

            for c in self.priv[2][i].coeff:
                s += str(c)+"|"
            s = s[:-1]
            s += '\n'
            s += "~"
        s = s[:-2]
        s += "\n\n-----END NTRU PRIVATE KEY BLOCK-----"
        if printk:
            print(s)
        return s

    def import_priv(self, s):
        self.name = ""
        self.email = ""
        cursor = 0
        while s[cursor] != '\n':
            cursor += 1
        cursor += 1
        while s[cursor] != '<':
            self.name += s[cursor]
            cursor += 1
        cursor += 1
        while s[cursor] != '>':
            self.email += s[cursor]
            cursor += 1
        for i in range(2):
            while (s[cursor] != '\n'):
                cursor += 1
            cursor += 1
        f = []
        fp = []
        h = []
        sn = ""
        while s[cursor] != '\n':
            if s[cursor] == '|':
                f.append(int(sn))
                sn = ""
            else:
                sn += s[cursor]
            cursor += 1
        f.append(int(sn))
        cursor += 1
        sn = ""
        while s[cursor] != '\n':
            if s[cursor] == '|':
                fp.append(int(sn))
                sn = ""
            else:
                sn += s[cursor]
            cursor += 1
        fp.append(int(sn))
        cursor += 1
        sn = ""
        while s[cursor] != '\n':
            if s[cursor] == '|':
                h.append(int(sn))
                sn = ""
            else:
                sn += s[cursor]
            cursor += 1
        h.append(int(sn))
        sn = ""
        self.N = len(f)
        F = pn.Polynome(N=self.N)
        F.coeff = np.array(f)
        Fp = pn.Polynome(N=self.N)
        Fp.coeff = np.array(fp)
        H = pn.Polynome(N=self.N)
        H.coeff = np.array(h)

        self.B = 1

        self.priv = ([F], [Fp], [H])



if __name__ == "__main__":
    import time
    t = time.time()
    N = 32
    k = KeyPair(N, 3*N//4, N//4, 5, 24, gen=True, name="Paul Martin", email="pmartin@email.fr")
    print(f"Time to calculate key : {int((time.time()-t)*100)/100}s")
    s = k.export(True)
    k2 = KeyPair()
    k2.import_pub(s)
    k2.export(True)
