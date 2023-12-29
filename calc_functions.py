from math import exp, log
from common import Point


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
    top = p + g0 + g1 * t + g2 * t**2
    bot = pref + g0 + g1 * t + g2 * t**2
    try:
        second_part = (top / bot) ** h
    except ZeroDivisionError:
        raise ValueError("Barometric part equals zero")
    return result * second_part


# point = [v1, v2], static_params = [y1, y2]
def refutas(point: list[float], static_params: list[float]) -> float:
    a = [14.534 * log(log(point[i] + 0.8)) + 10.975 for i in range(2)]
    a12 = sum([static_params[i] * a[i] for i in range(2)])
    return exp(exp((a12 - 10.975) / 14.534) - 0.8)


def arrhenius(point: list[Point]) -> float:
    try:
        return exp(sum([i.x * log(i.v) for i in point]))
    except TypeError:
        pass


def kendall_monroe(point: list[Point]) -> float:
    return (sum(i.x * i.v ** (1 / 3) for i in point)) ** 3


# params = d, point = [mu1, mu2, x1, x2]
def grunberg_nissan(params: list[float], point: list[Point]):
    d = [params[:2], params[2:]]
    visc_sum = sum([i.x * i.mu for i in point])
    param_sum = sum(
        [
            sum([point[i].x * point[j].mu * d[i][j] for j in range(len(point))])
            for i in range(len(point))
        ]
    )
    return visc_sum + param_sum


# params = alpha, point = [mu1, mu2, x1, x2]
def lederer_roegiers(params: list[float], point: list[Point]) -> float:
    return exp(
        (point[0].x * log(point[0].mu) + params[0] * point[1].x * log(point[1].mu))
        / (point[0].x + params[0] * point[1].x)
    )


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


def rho_full(
    params: list[float], static_params: list[float], point: list[float]
) -> float:
    e, beta = params
    rho0, p0, t0 = static_params
    p, t = point
    return (rho0 / (1 + beta * (t - t0))) / (1 - (p - p0) / e)


def rho_paredes_t(params: list[float], point: float) -> float:
    a = params
    return a[0] + a[1] * point + a[2] * point**2


def rho_paredes(
    params: list[float], static_params: list[float], point: list[float]
) -> float:
    c, b0, b1, b2 = params
    a0, a1, a2, pref = static_params
    p, t = point
    top = rho_paredes_t([a0, a1, a2], t)
    bottom_common = rho_paredes_t([b0, b1, b2], t)
    bottom_log = (bottom_common + p) / (bottom_common + pref)
    try:
        log(bottom_log)
    except ValueError:
        pass
    return top / (1 - c * log(bottom_log))
