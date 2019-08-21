"""
Byunghyun Ban
halfbottle@sangsang.farm
"""

import numpy as np
import math
import random
import time

dt = 5e-4
MSE = 1e-23

Vesicles = []


class Network():
    def __init__(self, filename, out_filename="result.csv"):
        self.vesicles = []
        self.nameidx = {}
        self.ks = {}
        self.flow_velocity = 0 # cm/s
        self._read_system(filename)
        self.Xs = np.zeros(len(self.vesicles))
        self.record = []
        self.out_filename = out_filename

        for i, el in enumerate(self.vesicles):
            self.Xs[i] = el.value

    def __repr__(self):
        return "system_simulator Network object with " + str(len(self.vesicles)) + " nodes"

    def _read_system(self, filename):
        network_file = open(filename, 'r')
        # reaction equation reading
        count = 0
        for line in network_file:
            if line[0] == "!":
                continue
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
                    self.vesicles.append(Vesicle(name, terms, initial_value, ion_molecular_weight, is_ion))
                else:
                    self.vesicles[self.nameidx[name]].terms += terms
            elif "VELOCITY" in line:
                self.flow_velocity = float(line.split("=")[1].strip())
            else:
                continue

    def _calc_gradient(self):
        gradient = np.zeros_like(self.Xs)
        for i, el in enumerate(self.vesicles):
            gradient[i] = el.calc_gradient(self.Xs, self.nameidx)
        return gradient

    def synchronous_update(self):
        gradient = self._calc_gradient()
        #self.record.append(np.copy(self.Xs))
        self.Xs += gradient
        for i in range(len(self.Xs)):
            if self.Xs[i] < 1e-100:
                self.Xs[i] = 0

    def converge(self, coef_multiply = False, save_logs = True):
        time_start = time.time()
        previous = np.copy(self.Xs)
        brk = False
        count = 0
        global dt

        if save_logs:
            res = open(self.out_filename, 'w')
            name = []
            for el in self.nameidx:
                name.append(el)
            res.write(", ".join(name) + "\n")

        while not brk:
            self.synchronous_update()
            for el in self.Xs:
                if math.isnan(el) or math.isinf(el):
                    self.Xs = previous
                    brk = True
                    continue

            error = self.Xs - previous
            mse = np.sum(error * error)

            if count % 100 == 0 and save_logs:
                res.write(str(self.Xs[0]))
                for i in range(len(self.Xs) - 1):
                    res.write(", " + str(self.Xs[i+1]))
                res.write("\n")

            if count % 50000 == 0 :
                print(str(count * dt / (1 + self.flow_velocity)) + " seconds past")
                print(str(count) + "th iteration\n>> current mse is : " + str(mse) + "\n")
                self.record.append(np.copy(self.Xs))
                self.show_result()
                if coef_multiply and mse < 1e-20:
                    dt *= 1.5

            count += 1
            if mse < MSE or math.isnan(mse):
                print("******** Simulation Done ********")
                self.show_result()
                print("The simulation took " + str(time.time() - time_start) + " seconds.")
                print("MSE is " + str(MSE))
                a = open(self.out_filename[:-4] + "txt")
                a.write("******** Simulation Done ********")
                self.write_result(a)
                a.close()
                break
            previous = np.copy(self.Xs)

    def show_result(self):
        for name in self.nameidx:
            if not self.vesicles[self.nameidx[name]].is_ion:
                    print("%15s@ : "%name + "  " + str(self.Xs[self.nameidx[name]]) + " mol/L")
            else:
                print("%15s  : "%name + "  " + str(self.Xs[self.nameidx[name]]) + " mol/L")
        print("pH is\t: " + str(-math.log10(self.Xs[self.nameidx["[H+]"]])))
        print("TDS is\t: " + str(self.calc_ppm()) + " mg/L(ppm)\n")

    def write_result(self, file):
        for name in self.nameidx:
            if not self.vesicles[self.nameidx[name]].is_ion:
                    file.write("%15s@ : "%name + "  " + str(self.Xs[self.nameidx[name]]) + " mol/L")
            else:
                file.write("%15s  : "%name + "  " + str(self.Xs[self.nameidx[name]]) + " mol/L")
        file.write("pH is\t: " + str(-math.log10(self.Xs[self.nameidx["[H+]"]])))
        file.write("TDS is\t: " + str(self.calc_ppm()) + " mg/L(ppm)\n")

    def calc_ppm(self):
        total_mass = 0
        for el in self.vesicles:
            if el.is_ion:
                total_mass += el.molecular_weight * el.value
        return total_mass * 1000

    def export_cytoscape(self, filename="cytoscape_form_export.csv"):
        if filename[-4:] != ".csv":
            filename += ".csv"
        res = open(filename, 'w')
        res.write("Target,Source,Interaction")
        target_source_interaction = []
        for vesicle in self.vesicles:
            for term in vesicle.terms:
                for elements in term.elements:
                    if term.coefficient == 0:
                        continue
                    elif term.coefficient > 0:
                        interaction = "1"
                    else:
                        interaction = "-1"

                    if (vesicle.name[1:-1], elements[1:-1], interaction) in target_source_interaction:
                        continue
                    res.write(",".join((vesicle.name[1:-1], elements[1:-1], interaction)))
                    res.write("\n")
                    target_source_interaction.append((vesicle.name[1:-1], elements[1:-1], interaction))
        res.close()

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
    def __init__(self, name, terms, initial_value=0, molecular_weight=0, is_ion=True):
        self.name = name
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
        if self.value < 1e-100:
            self.value = 0
        return dX
