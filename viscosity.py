import matplotlib.pyplot as plt
from numpy import arange, array, zeros
from scipy.optimize import minimize
import matplotlib
from thermo import Chemical

from calc_functions import (
    DEF,
    HG,
    grunberg_nissan,
    lederer_roegiers,
    arrhenius,
    kendall_monroe,
    refutas,
)
from common import (
    Point,
    ViscData,
    aad,
    aad_new,
    create_frac_calc,
    create_frac_calc_2comp,
    normalize,
    mole_to_mass,
)
from data_csv import get_data
from Flashers_and_classes.classes import FlashCalculator
from density import Density
from material import Material


class VFT:
    def __init__(self):
        self.pref = 0
        self.params = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.maxiter = 2**20

    def setup_def(self, mu_pref: dict[float:float]) -> bool:
        visc_pref = list(mu_pref.values())
        temps = list(mu_pref.keys())

        # Выбираем начальное приближение на основе характера экспериментальных данных
        x0 = [0.000852102519255135, -67.82280822979249, 4.052812591576497]

        args = (DEF, temps, visc_pref)
        minimized = minimize(
            fun=aad_new,
            x0=array(x0),
            method="Nelder-Mead",
            args=args,
            tol=1e-5,
            options={"maxiter": self.maxiter, "maxfev": self.maxiter},
        )
        self.params[0], self.params[1], self.params[2] = minimized.x
        print(minimized.message)
        return minimized.success

    def setup_hg(self, pref: float, mu_data: dict[float : dict[float:float]]) -> bool:
        points = []
        exp_visc = []
        p_arr = list(mu_data.keys())
        for p in p_arr:
            t_arr = list(mu_data[p].keys())
            for t in t_arr:
                points.append([p, t])
                exp_visc.append(mu_data[p][t])

        x0 = [-1663, 10, -0.01, 7]

        args = (HG, points, exp_visc, self.params[:3] + [pref])
        minimized = minimize(
            fun=aad_new,
            x0=array(x0),
            method="Nelder-Mead",
            args=args,
            tol=1e-7,
            options={"maxiter": self.maxiter, "maxfev": self.maxiter},
        )
        self.params[3], self.params[4], self.params[5], self.params[6] = minimized.x
        return minimized.success

    def setup(self, pref: float, mu_data: dict[float : dict[float:float]]) -> bool:
        self.pref = pref
        def_flag = self.setup_def(mu_data[pref])
        if not def_flag:
            return False
        return self.setup_hg(pref, mu_data)

    def calc_visc_def(self, t: float) -> float:
        return DEF(self.params[:3], t)

    def calc_visc(self, t: float, p: float) -> float:
        return HG(self.params[3:], self.params[:3] + [self.pref], [p, t])


class LedererRoegiers:
    def __init__(
        self,
        calc_method: callable,
        vft_parameters: list[list[float]],
        frac_solver: FlashCalculator,
        x0: float,
    ):
        self.calc_model = calc_method
        self.frac_solver = frac_solver
        self.components_amount = len(vft_parameters)
        self.vft_solvers = [VFT() for i in range(self.components_amount)]
        for i in range(self.components_amount):
            self.vft_solvers[i].params = vft_parameters[i]
        self.params = [x0 for i in range(self.components_amount - 1)]

    def calc(self, t: float, p: float):
        visc_parts = [
            self.vft_solvers[i].calc_visc(t, p) for i in range(self.components_amount)
        ]
        moles = self.frac_solver.calculate_flash(t, p)
        last_visc = lederer_roegiers(
            self.params[0], normalize(moles[:2]), visc_parts[2:]
        )
        last_mole = sum(moles[:2])
        partitions = self.components_amount - 2
        for i in range(partitions):
            last_visc = lederer_roegiers(self.params[i + 1])


def main_visc_check():
    r22_chem = Chemical("r-22")
    r22_mat = Material(r22_chem)
    r115_chem = Chemical("Chloropentafluoroethane")
    r115_mat = Material(r115_chem)
    substance = "R22"
    data = ViscData(substance)
    pref = 0.06

    pref_visc = {i.t: i.mu for i in data.get_by_pressure(pref)}
    # oil_visc_d = [412.33, 136.68, 57.08, 28.17, 16.01, 9.83, 6.69, 4.73]
    # visc = {0.06: {t[i]: oil_visc_d[i] for i in range(len(t))}}
    temps = list(pref_visc.keys())
    viscs = [pref_visc[i] for i in temps]

    test = VFT()
    test.setup_def(pref_visc)
    print(
        f"D = {test.params[0]}\nE = {test.params[1]}\nF = {test.params[2]}\nG = {test.params[3:6]}\nH = {test.params[6]}"
    )

    font = {"family": "Times New Roman", "size": 20}
    matplotlib.rc("font", **font)
    fig, ax = plt.subplots()
    ax.scatter(temps, viscs, label="Экспериментальные данные")
    viscs_calc = [test.calc_visc(i, pref) for i in arange(temps[0], temps[-1] + 1)]

    mat_values = [r115_mat.calc_prop(i).mu for i in arange(temps[0], temps[-1] + 1)]
    viscs_temp = [test.calc_visc_def(i) for i in arange(temps[0], temps[-1] + 1)]
    ax.plot(
        arange(temps[0], temps[-1] + 1), viscs_calc, "r", label="Расчётные значения"
    )
    ax.plot(arange(temps[0], temps[-1] + 1), mat_values, "g", label="Thermo")
    ax.set_ylabel("$\it{Вязкость, Па*С}$")
    ax.set_xlabel("$\it{T, °C}$")
    plt.grid()
    plt.legend()
    plt.show()


def oil_coeff_calc():
    t_arr = [-25, 0, 25, 50, 75, 100, 125]
    t_calc = arange(t_arr[0], t_arr[-1] + 0.1, 0.1)
    rho_arr = [1733.38, 289.693, 64.8321, 19.9096, 7.79992, 4.19357, 3.09385]
    oil_calc = VFT()
    oil_calc.setup_def({t_arr[i]: rho_arr[i] for i in range(len(t_arr))})
    print(oil_calc.params)
    font = {"family": "Times New Roman", "size": 20}
    matplotlib.rc("font", **font)
    fig, ax = plt.subplots()
    ax.scatter(t_arr, rho_arr, label="Экспериментальные данные")
    ax.plot(
        t_calc,
        [oil_calc.calc_visc_def(i) for i in t_calc],
        "r",
        label="Расчётные значения",
    )
    ax.set_ylabel("$\it{Вязкость, сП}$")
    ax.set_xlabel("$\it{T, °C}$")
    plt.yscale("log")
    plt.grid()
    plt.legend()
    plt.show()


if __name__ == "__main__":
    oil_coeff_calc()
