from math import exp, log


# Params = [D, E, F], point = t
def DEF(params: list[float], point: float):
    d, e, f = params
    return d * exp(e / (point + 273.15 - f))


# Params = [G0, G1, G2, H], Static_Params = все для DEF и pref, point = [p, t]
def HG(params: list[float], static_params: list[float], point: list[float]):
    g0, g1, g2, h = params
    d, e, f, pref = static_params
    p, t = point
    result = DEF([d, e, f], t)
    top = (p + g0 + g1 * t + g2 * t ** 2)
    bot = (pref + g0 + g1 * t + g2 * t ** 2)
    second_part = (top / bot) ** h
    return result * second_part


# Params = d, Static_Params = x, point = mu
def grunberg_nissan(params: list[list[float]], static_params: list[float], point: list[float]):
    d = params
    x = static_params
    mu = point
    visc_sum = sum([x[i] * log(mu[i]) for i in range(2)])
    param_sum = sum([sum([x[i] * x[j] * d[i][j]] for j in range(2)) for i in range(2)])
    return visc_sum + param_sum
