from common import aad
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from os import listdir
from numpy import arange


class Density:
    def __init__(self, rho0, t0, p0):
        self.rho0 = rho0
        self.t0 = t0
        self.p0 = p0
        self.beta = 0
        self.e = 10

    def calc_rho1(self, t, beta=None):
        if not beta:
            beta = self.beta
        return self.rho0 / (1 + beta * (t - self.t0))

    def calc_rho2(self, p, e=None):
        if not e:
            e = self.e
        return self.rho0 / (1 - (p - self.p0) / e)

    def calc_rho_full(self, p, t, e=None, beta=None):
        if not e:
            e = self.e
        if not beta:
            beta = self.beta
        return (self.rho0 / (1 + beta * (t - self.t0))) / (1 - (p - self.p0) / e)

    def aad_rho1(self, beta, temp, rho):
        rho_calc = [self.calc_rho1(i, beta) for i in temp]
        return aad(rho, rho_calc)

    def aad_rho2(self, e, p, rho):
        rho_calc = [self.calc_rho2(i, e) for i in p]
        return aad(rho, rho_calc)

    def aad_rho_full(self, params, rho):
        e = params[0]
        beta = params[1]
        rho_calc = []
        rho_exp = []
        for p_cur in list(rho.keys()):
            for t_cur in list(rho[p_cur].keys()):
                rho_calc.append(self.calc_rho_full(p_cur, t_cur, e, beta))
                rho_exp.append(rho[p_cur][t_cur])
        test = aad(rho_calc, rho_exp)
        return test

    def setup_rho1(self, temps, rho):
        args = (temps, rho)
        test = minimize(fun=self.aad_rho1, x0=0, args=args)
        self.beta = test.x[0]

    def setup_rho2(self, p, rho):
        args = (p, rho)
        test = minimize(fun=self.aad_rho2, x0=3, args=args)
        self.e = test.x[0]

    def setup_rho_full(self, rho):
        args = (rho)
        test = minimize(fun=self.aad_rho_full, method='Nelder-Mead', x0=[3, 0], args=args)
        self.e = test.x[0]
        self.beta = test.x[1]


if __name__ == "__main__":
    rho = {}
    substance = "R22"
    for file in listdir(substance):
        with open(substance + "\\" + file, "r") as f:
            text = f.read()
            p = float(file[:-4])
            rho[p] = {}
            for line in text.split('\n'):
                cur_line = line.split(' ')
                rho[p][float(cur_line[0])] = float(cur_line[2])
    p = list(rho.keys())
    graph_numb = [0, 0]
    p0 = 0.05
    temps = list(rho[p0].keys())
    t0 = -80
    rho_p_dynamic = [rho[i][t0] for i in p]
    rho_t = list(rho[p0].values())
    rho_p = [rho[i][t0] for i in p]
    rho0 = rho[p0][t0]
    test = Density(rho0, t0, p0)
    test.setup_rho_full(rho)
    fig, ax = plt.subplots(2, 4)

    t_exp = list(rho[p0].keys())
    t_calc = arange(t_exp[0], t_exp[-1], 0.5)
    for i in range(len(p)):
        p_cur = p[i]
        graphy = (i + 1) // 4
        graphx = (i + 1) % 4
        exp_temps = list(rho[p_cur].keys())
        t_arr = arange(exp_temps[0], exp_temps[-1], 0.5)
        ax[graphy][graphx].plot(t_arr, [test.calc_rho_full(p_cur, j) for j in t_arr], "r", label="Calculated")
        ax[graphy][graphx].plot(exp_temps, list(rho[p_cur].values()), "b", label="Expected")
        ax[graphy][graphx].set_title("P = %f" % p_cur)
    print("Beta = %f" % test.beta)
    print("E = %f" % test.e)
    print(test.aad_rho2(test.e, p, rho_p))
    plt.legend()
    plt.show()
