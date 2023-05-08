from common import aad, ViscData, aad_new
from scipy.optimize import minimize, Bounds
import matplotlib.pyplot as plt
from numpy import arange
from calc_functions import DEF, HG, grunberg_nissan


class VFT:
    def __init__(self, mu_data: dict[float:dict[float:float]], pref: float):
        self.mu_data = mu_data
        self.pref = pref
        self.params = [0, 0, 0, 0, 0, 0, 0]

    def setup_def(self):
        visc_pref = list(self.mu_data[pref].values())
        temps = list(self.mu_data[pref].keys())

        # Выбираем начальное приближение на основе характера экспериментальных данных
        x0 = [0.010124325475473392, 498.84631096931736, 372.92721289829717]
        if visc_pref[0] < visc_pref[-1]:
            x0 = [8.956959670631126e-05, -1865.4191615857121, 1178.868480983142]

        args = (DEF, temps, visc_pref)
        test = minimize(fun=aad_new, x0=x0, method='Nelder-Mead', args=args, tol=1e-5)
        print(test.message)
        self.params[0], self.params[1], self.params[2] = test.x

    def setup_hg(self):
        points = []
        exp_visc = []
        p_arr = list(self.mu_data.keys())
        for p in p_arr:
            t_arr = list(self.mu_data[p].keys())
            for t in t_arr:
                points.append([p, t])
                exp_visc.append(self.mu_data[p][t])

        x0 = [-1663, 10, -0.01, 7]

        args = (HG, points, exp_visc, self.params[:3] + [self.pref])
        test = minimize(fun=aad_new, x0=x0, method='Nelder-Mead', args=args, tol=1e-7)
        self.params[3], self.params[4], self.params[5], self.params[6] = test.x

    def calc_visc_def(self, t):
        return DEF(self.params[:3], t)

    def calcVisc(self, t, p):
        return HG(self.params[3:], self.params[:3] + [self.pref], [p, t])


class MuMix:
    def __init__(self, ):
        ...


if __name__ == "__main__":
    substance = "R115"
    data = ViscData(substance)
    pref = 0.06

    visc = data.visc_d
    temps = list(visc[pref].keys())
    viscs = [i[2] for i in data.get_by_pressure(2, pref)]

    test = VFT(visc, pref)
    test.setup_def()
    test.setup_hg()
    print(f'D = {test.params[0]}\nE = {test.params[1]}\nF = {test.params[2]}\n')
    # print(f'AAD = {aadDEF(test.params[:3], temps, viscs)}')

    fig, ax = plt.subplots()
    ax.scatter(temps, viscs, label="Экспериментальные данные")
    viscs_calc = [test.calcVisc(i, pref) for i in arange(temps[0], temps[-1] + 1)]
    viscs_temp = [test.calc_visc_def(i) for i in arange(temps[0], temps[-1] + 1)]
    ax.plot(arange(temps[0], temps[-1] + 1), viscs_calc, "r",
            label="Расчётные значения")
    ax.set_ylabel("Mu, cP")
    ax.set_xlabel("T, Гр. Цельсия")
    ax.set_title("mu = f(t)")
    plt.grid()
    plt.legend()
    plt.show()
