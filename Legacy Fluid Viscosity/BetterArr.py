from numpy import multiply, divide, add, subtract, negative, exp, empty, log
from array import array


class BetterArr(array):
    def __new__(cls, *args, **kwargs):
        return super(BetterArr, cls).__new__(cls, 'f', *args, **kwargs)

    def __mul__(self, other):
        return multiply(self, other)

    def __truediv__(self, other):
        return divide(self, other)

    def __add__(self, other):
        return add(self, other)

    def __sub__(self, other):
        return subtract(self, other)

    def __neg__(self):
        return negative(self)
