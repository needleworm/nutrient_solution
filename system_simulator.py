import numpy as np

dt = 0.1

Vesicles = []


class Network():
    def __init__(self, filename):
        self.vesicles = []
        self.nameidx = {}
        self.ks = {}
        self._read_system(filename)
        self.Xs = np.zeros(len(self.vesicles))
        self.record = []

        for i, el in enumerate(self.vesicles):
            self.Xs[i] = el.value

    def __repr__(self):
        return "system_simulator Network object with " + str(len(self.vesicles)) + " nodes"

    def _read_system(self, filename):
        network_file = open(filename, 'r')
        network_file.readline() # remove header

        # reaction rate constant reading
        for line in network_file:
            if line[0] == 'k':
                k, value = line.strip().split(" = ")
                self.ks[k] = float(value.strip())
            else:
                break

        # reaction equation reading
        count = 0
        for line in network_file:
            if line[0] != "&":
                break

            initial_value = 0

            if "$" in line:
                line, value = line.split("$")
                try:
                    initial_value = float(value.strip())
                except:
                    print (line + "  has wrong initial_mol indication")
                    exit(1)

            splt = line.split("*")
            name = splt[0][1:].strip()
            self.nameidx[name] = count
            count += 1

            splt = splt[1:]
            terms = []
            for el in splt:
                if "k" not in el:
                    print(line + " has wrong term " + el)
                    exit(1)

                tokens = el.strip().split("[")
                coef = tokens[0]
                if "-" in coef:
                    coefficient = -1 * self.ks[coef[1:]]
                else:
                    coefficient = self.ks[coef]

                elements = []
                for el in tokens[1:]:
                    elements.append("[" + el)

                terms.append(Terms(coefficient, elements))

            self.vesicles.append(Vesicle(terms, initial_value))

    def _calc_gradient(self):
        gradient = np.zeros_like(self.Xs)
        for i, el in enumerate(self.vesicles):
            gradient[i] = el.calc_gradient(self.Xs, self.nameidx)
        return gradient

    def synchronous_update(self):
        gradient = self._calc_gradient()
        self.record.append(np.copy(self.Xs))
        self.Xs += gradient

class Terms():
    def __init__(self, coefficient, elements):
        self.coefficient = coefficient
        self.elements = elements # (idx_1, idx_2, idx_3 ... idx_n)

    def __repr__(self):
        return str(self.coefficient) + " " + str(self.elements)

    def term(self, Xs, nameidx):
        retval = self.coefficient
        for el in self.elements:
            retval *= Xs[nameidx[el]]
        return retval


class Vesicle():
    def __init__(self, terms, initial_value=0):
        self.value = initial_value
        self.terms = terms

    def __repr__(self):
        return str(self.terms)

    def calc_gradient(self, Xs, nameidx):
        dX_dt = 0
        for term in self.terms:
            dX_dt += term.term(Xs, nameidx)

        dX = dX_dt * dt
        self.value += dX
        return dX
