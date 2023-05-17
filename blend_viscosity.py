from math import log
from dataclasses import dataclass
from common import aad
from scipy.optimize import minimize
from os import listdir
import matplotlib.pyplot as plt
from numpy import arange

class BlendVisc:
    def __init__(self, x: list[float], mu: list[float]):
        pass