import Polynome as pn
import numpy as np


class KeyPair:
    def __init__(self,
                 N=256,
                 p=2,
                 q=127,
                 d=128,
                 Bk=8,
                 gen=False,
                 name="User Name",
                 email="user@example.com"):
        """
        Create a key with the parameter passed to the constructor
        """
        if gen:
            # Find f that is inversible modulo q
            while True:
                try:
                    f = pn.T(d, N)
                    invf = f.inv(q)
                    if pn.NormF(f) < Bk:
                        break
                except Exception:
                    pass
            # Find g that is inversible modulo p
            while True:
                try:
                    g = pn.T(d, N)
                    g.inv(p)
                    if pn.NormF(g) < Bk:
                        break
                except Exception:
                    pass

            # h = g/fp (mod q)
            h = pn.modXnp1(g * invf, g.N)
            for i in range(len(h)):
                h.coeff[i] = pn.xgcd(h.coeff[i], q)[0]

            self.pub = h
            self.priv = (f*p, g)

        self.N = N
        self.Bk = Bk
        self.q = q
        self.d = d
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
        for i in range(self.B+1):
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
            s += "~\n"
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

        self.priv = [[], [], []]

        self.B = 0
        while True:
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
            self.priv[0].append(F)
            self.priv[1].append(Fp)
            self.priv[2].append(H)

            cursor += 1
            if s[cursor] != '~':
                break
            cursor += 2
            self.B += 1


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
