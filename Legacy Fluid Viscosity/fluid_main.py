from BetterArr import BetterArr
from numpy import empty, exp, log

def eps_calc(kb, tc):
    return kb * tc / 1.2593


def sigma_calc(vc, na):
    return (0.302 * vc / na) ** (1.0 / 3.0)


def sigmaij_calc(sigma_ii, sigma_jj):
    return ((sigma_ii ** 3 + sigma_jj ** 3) / 2) ** (1 / 3)


def sigma_setup(vc, na):
    sigmamain = empty(len(vc))
    for i in range(len(vc)):
        sigmamain[i] = sigma_calc(vc[i], na)
    returnable = empty([len(vc), len(vc)])
    for i in range(len(vc)):
        returnable[i][i] = sigmamain[i]
        for j in range(i + 1, len(vc)):
            returnable[i][j] = returnable[j][i] = sigmaij_calc(sigmamain[i], sigmamain[j])
    return returnable


def eps_setup(sigma, kb, tc):
    eps = empty([len(sigma), len(sigma)])
    for i in range(len(sigma)):
        eps[i][i] = eps_calc(kb, tc[i])
    for i in range(len(sigma)):
        for j in range(i + 1, len(sigma)):
            eps[i][j] = eps[j][i] = (eps[i][i] * sigma[i][i] ** 3 + eps[j][j] * sigma[j][j] ** 3) / 2.0 / sigma[i][
                j] ** 3
    return eps


# вынести tc и vc в внешний вектор
def viscosity(x, mw_i, nu_i, t):
    kb = 1.380649 * 1e-23
    na = 6.022 * 1e23
    tc = [516.25, 540.61]
    vc = [182.49 * 1e-6, 428.64 * 1e-6]
    sigma = sigma_setup(vc, na)
    eps = eps_setup(sigma, kb, tc)
    sigma_x = 0.0
    eps_x = 0.0
    nu = 0.0
    mwx = 0.0
    for i in range(0, 2):
        nu += nu_i[i] * x[i]
        mwx += mw_i[i] * x[i]
        for j in range(0, 2):
            sigma_x += x[i] * x[j] * sigma[j][i] ** 3
            eps_x += x[i] * x[j] * eps[j][i] * sigma[j][i] ** 3
    sigma_x = sigma_x ** (1.0 / 3.0)
    eps_x /= sigma_x ** 3
    t_ast = kb * t / eps_x
    rho_ast = na * (sigma_x ** 3) / nu
    # не меняется
    coeff = [1.06036, 0.15610, 0.19300, 0.47635, 1.03587, 1.52996, 1.76474, 3.89411]
    # можно подогнать
    b = [0.062692, 4.095577, -8.743269 * 1e-6, 11.124920, 2.542477 * 1e-6, 14.863984]

    collision_integral = coeff[0] / t_ast ** coeff[1] + coeff[2] / exp(coeff[3] * t_ast) + coeff[4] / exp(
        coeff[5] * t_ast) + coeff[6] / exp(coeff[7] * t_ast)
    ac = 0.95
    visc_0_ast = 0.17629 * t_ast ** (1.0 / 2.0) * ac / collision_integral
    delta_visc_ast = b[0] * (exp(b[1] * rho_ast) - 1) + b[2] * (exp(b[3] * rho_ast) - 1) + b[4] * t_ast ** (-2) * (
                exp(b[5] * rho_ast) - 1)
    visc_ast = visc_0_ast + delta_visc_ast
    visc = visc_ast * (mwx * eps_x / na) ** 0.5 / (sigma_x ** 2)
    return visc


def LJ_fluid():
    vexp = [
        BetterArr([0.410, 0.419, 0.450, 0.494, 0.556, 0.644, 0.757, 0.931, 1.194]),
        BetterArr([0.509, 0.518, 0.547, 0.597, 0.668, 0.760, 0.885, 1.052, 1.344]),
        BetterArr([0.604, 0.614, 0.648, 0.707, 0.787, 0.888, 1.023, 1.200, 1.494]),
        BetterArr([0.706, 0.718, 0.757, 0.824, 0.913, 1.022, 1.168, 1.354, 1.647]),
        BetterArr([0.817, 0.832, 0.875, 0.951, 1.047, 1.162, 1.316, 1.509, 1.797]),
        BetterArr([0.937, 0.953, 1.003, 1.087, 1.191, 1.312, 1.473, 1.674, 1.955])]
    T = 293.15
    P = BetterArr([0.1, 20, 40, 60, 80, 100])
    Mw = BetterArr([46.07, 100.2]) / 1000
    constC = 0.0894

    # настроечные параметры
    B0 = BetterArr([1e2, 1e2])
    B1 = BetterArr([0.0009, 0.0009])
    alpha_n = BetterArr([0.002, 0.004])
    nu0 = BetterArr([1.04958e-4, 4.76409e-4])
    nuP0 = nu0 * exp(-alpha_n * T)

    B = B0 * exp(-B1 * T)
    vvv = empty((len(P), len(vexp)))
    for j in range(0, len(P)):
        nu = nuP0 * (1 - constC * log(1 + P[j] * (1 / B)))
        for i in range(0, len(vexp)):
            x1 = (i - 1) * 0.125
            x2 = 1 - x1
            vvv[i][j] = viscosity([x1, x2], Mw, nu, T)
    return vvv


if __name__ == "__main__":
    output = LJ_fluid()
    print(output)
