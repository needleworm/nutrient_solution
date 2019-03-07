import numpy as np
import math

dt = float("1e-4")

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

        # reaction equation reading
        count = 0
        for line in network_file:
            if line[0] == 'k':
                k, value = line.strip().split("=")
                self.ks[k.strip()] = float(value.strip())
            elif line[0] == "&":

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

                if name not in self.nameidx:
                    self.nameidx[name] = count
                    count += 1
                    self.vesicles.append(Vesicle(terms, initial_value))
                else:
                    self.vesicles[self.nameidx[name]].terms += terms
            else:
                continue

    def _calc_gradient(self):
        gradient = np.zeros_like(self.Xs)
        for i, el in enumerate(self.vesicles):
            gdnt = el.calc_gradient(self.Xs, self.nameidx)
            gradient[i] = gdnt
        return gradient

    def synchronous_update(self):
        gradient = self._calc_gradient()
        self.record.append(np.copy(self.Xs))
        self.Xs += gradient
        for i in range(len(self.Xs)):
            if self.Xs[i] < float("1e-32") :
                self.Xs[i] = 0

    def converge(self):
        previous = np.copy(self.Xs)
        brk = False
        count = 0
        while not brk:
            self.synchronous_update()
            for el in self.Xs:
                if math.isnan(el) or math.isinf(el):
                    self.Xs = previous
                    brk = True
                    continue
            error = self.Xs - previous
            mse = np.sum(error * error)
            if count % 1000 == 0 :
                print(str(count) + "th iteration\n>> current mse is : " + str(mse))
            count += 1
            if mse < float("1e-21"):
                break
            previous = np.copy(self.Xs)

    def show_result(self):
        for name in self.nameidx:
            print(name + "\t: " + str(self.Xs[self.nameidx[name]]))
        print("pH is\t: " + str(-math.log10(float("1e-7") + self.Xs[self.nameidx["[H+]"]])))

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
        #if abs(dX) < float("1e-16"):
        #    dX = 0
        self.value += dX
        return dX
