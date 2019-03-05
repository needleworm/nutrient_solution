import numpy as np

dt = 0.1

Vesicles = []


class Network():
    def __init__(self, initial_Xs):
        self.vesicles = []


class Terms():
    def __init__(self, coefficient, elements):
        self.coefficient = coefficient
        self.elements = elements # (idx_1, idx_2, idx_3 ... idx_n)

    def term(self, Xs):
        retval = self.coefficient
        for el in self.elements:
            retval *= Xs[el]
        return retval


class Vesicle():
    def __init__(self, initial_value, idx_X, terms):
        self.value = initial_value
        self.idx = idx_X
        self.terms = terms

    def calc_gradient(self, Xs):
        dX_dt = 0
        for term in self.terms:
            dX_dt += term.term(Xs)

        dX = dX_dt * dt
        self.value += dX
        return dX
