import numpy as np
import math
import random

dt = 1e-5
#Avogadro_constant = 6.02214129e23
Avogadro_constant = 1

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
                ion_molecular_weight = 0
                is_ion = False
                if "#" in line:
                    line = line.split("#")[0].strip()
                    is_ion = True

                if "@" in line:
                    line, value = line.split("@")
                    ion_molecular_weight = float(value.strip())

                if "$" in line:
                    line, value = line.split("$")
                    try:
                        initial_value = float(value.strip()) * Avogadro_constant
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
                    self.vesicles.append(Vesicle(terms, initial_value, ion_molecular_weight, is_ion))
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
        #self.record.append(np.copy(self.Xs))
        self.Xs += gradient
        for i in range(len(self.Xs)):
            if self.Xs[i] < 1e-300 :
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
            if count % 100000 == 0 :
                print(str(count * dt) + " seconds past")
                print(str(count) + "th iteration\n>> current mse is : " + str(mse) + "\n")
                self.record.append(np.copy(self.Xs))
                self.show_result()
            count += 1
            if mse < 1e-19 or math.isnan(mse):
                break
            previous = np.copy(self.Xs)

    def show_result(self):
        for name in self.nameidx:
            if not self.vesicles[self.nameidx[name]].is_ion:
                    print("%15s@ : "%name + "  " + str(self.Xs[self.nameidx[name]]/Avogadro_constant) + " mol/L")
            else:
                print("%15s  : "%name + "  " + str(self.Xs[self.nameidx[name]]/Avogadro_constant) + " mol/L")
        print("pH is\t: " + str(-math.log10(self.Xs[self.nameidx["[H+]"]])))
        print("TDS is\t: " + str(self.calc_ppm()) + " mg/L(ppm)\n")

    def calc_ppm(self):
        total_mass = 0
        for el in self.vesicles:
            if not el.is_ion:
                continue
            total_mass += el.molecular_weight * el.value
        return total_mass *1000

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
    def __init__(self, terms, initial_value=0, molecular_weight=0, is_ion=True):
        self.value = initial_value
        self.terms = terms
        self.molecular_weight = molecular_weight
        self.is_ion = is_ion

    def __repr__(self):
        return str(self.terms)

    def calc_gradient(self, Xs, nameidx):
        dX_dt = 0
        for term in self.terms:
            dX_dt += term.term(Xs, nameidx)

        dX = dX_dt * dt
        self.value += dX
        return dX
