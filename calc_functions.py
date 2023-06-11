from math import exp, log


# Params = [D, E, F], point = t
def DEF(params: list[float], point: float) -> float:
    d, e, f = params
    return d * exp(e / (point + 273.15 - f))


# Params = [G0, G1, G2, H], Static_Params = все для DEF и pref, point = [p, t]
def HG(params: list[float], static_params: list[float], point: list[float]) -> float:
    g0, g1, g2, h = params
    d, e, f, pref = static_params
    p, t = point
    result = DEF([d, e, f], t)
    top = (p + g0 + g1 * t + g2 * t ** 2)
    bot = (pref + g0 + g1 * t + g2 * t ** 2)
    second_part = (top / bot) ** h
    return result * second_part


# point = [v1, v2], static_params = [y1, y2]
def refutas(point: list[float], static_params: list[float]) -> float:
    a = [14.534 * log(log(point[i] + 0.8)) + 10.975 for i in range(2)]
    a12 = sum([static_params[i] * a[i] for i in range(2)])
    return exp(exp((a12 - 10.975) / 14.534) - 0.8)


# point = [mu1, mu2], static_params: x1, x2
def arrhenius(point: list[float], static_params: list[float]) -> float:
    mu1, mu2 = point
    x1, x2 = static_params
    return exp(x1 * log(mu1) + x2 * log(mu2))


# point = [v1, v2], static_params = [x1, x2]
def kendall_monroe(point: list[float], static_params: list[float]) -> float:
    v1, v2 = point
    x1, x2 = static_params
    return (x1 * v1 ** (1 / 3) + x2 * v2 ** (1 / 3)) ** 3


# params = d, static_params = x, point = [x1, x2, mu1, mu2]
def grunberg_nissan(params: list[float], point: list[float]):
    d = [params[:2], params[2:]]
    x = point[2:]
    mu = point[:2]
    visc_sum = sum([x[i] * log(mu[i]) for i in range(2)])
    param_sum = sum([sum([x[i] * x[j] * d[i][j]] for j in range(2)) for i in range(2)])
    return visc_sum + param_sum


def lederer_roegiers(params: float, static_params: list[float], point: list[float]) -> float:
    alpha = params
    x1, x2 = static_params
    mu1, mu2 = point
    return exp((x1 * log(mu1) + alpha * x2 * log(mu2)) / (x1 + alpha * x2))


def rho_t(params: float, static_params: list[float], point: float) -> float:
    e = params
    rho0, t0 = static_params
    t = point
    return rho0 / (1 + e * (t - t0))


def rho_p(params: float, static_params: list[float], point: float) -> float:
    beta = params
    rho0, p0 = static_params
    p = point
    return rho0 / (1 - (p - p0) / beta)


def rho_full(params: list[float], static_params: list[float], point: list[float]) -> float:
    e, beta = params
    rho0, p0, t0 = static_params
    p, t = point
    return (rho0 / (1 + beta * (t - t0))) / (1 - (p - p0) / e)
