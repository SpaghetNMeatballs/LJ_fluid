import matplotlib.pyplot as plt
from numpy import arange, array, zeros
from scipy.optimize import minimize
import matplotlib
from calc_functions import DEF, HG, grunberg_nissan, lederer_roegiers, arrhenius, kendall_monroe, refutas
from common import ViscData, aad, aad_new, create_frac_calc, normalise, mole_to_mass
from data_csv.cavestri_loader import get_data
from Flashers_and_classes.classes import FlashCalculator


class VFT:
    def __init__(self):
        self.pref = 0
        self.params = [0., 0., 0., 0., 0., 0., 0.]
        self.maxiter = 2 ** 15

    def setup_def(self, pref: float, mu_data: dict[float:dict[float:float]]) -> bool:
        visc_pref = list(mu_data[pref].values())
        temps = list(mu_data[pref].keys())

        # Выбираем начальное приближение на основе характера экспериментальных данных
        x0 = [0, 0, 0]

        args = (DEF, temps, visc_pref)
        minimized = minimize(fun=aad_new, x0=array(x0), method='Nelder-Mead', args=args, tol=1e-5,
                             options={'maxiter': self.maxiter, 'maxfev': self.maxiter})
        self.params[0], self.params[1], self.params[2] = minimized.x
        print(minimized.message)
        return minimized.success

    def setup_hg(self, pref: float, mu_data: dict[float:dict[float:float]]) -> bool:

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
        minimized = minimize(fun=aad_new, x0=array(x0), method='Nelder-Mead', args=args, tol=1e-7,
                             options={'maxiter': self.maxiter, 'maxfev': self.maxiter})
        self.params[3], self.params[4], self.params[5], self.params[6] = minimized.x
        return minimized.success

    def setup(self, pref: float, mu_data: dict[float:dict[float:float]]) -> bool:
        self.pref = pref
        def_flag = self.setup_def(pref, mu_data)
        if not def_flag:
            return False
        return self.setup_hg(pref, mu_data)

    def calc_visc_def(self, t: float) -> float:
        return DEF(self.params[:3], t)

    def calc_visc(self, t: float, p: float) -> float:
        return HG(self.params[3:], self.params[:3] + [self.pref], [p, t])


class LedererRoegiers:
    def __init__(self, calc_method: callable, vft_parameters: list[list[float]], frac_solver: FlashCalculator,
                 x0: float):
        self.calc_model = calc_method
        self.frac_solver = frac_solver
        self.components_amount = len(vft_parameters)
        self.vft_solvers = [VFT() for i in range(self.components_amount)]
        for i in range(self.components_amount):
            self.vft_solvers[i].params = vft_parameters[i]
        self.params = [x0 for i in range(self.components_amount - 1)]

    def calc(self, t: float, p: float):
        visc_parts = [self.vft_solvers[i].calc_visc(t, p) for i in range(self.components_amount)]
        moles = self.frac_solver.calculate_flash(t, p)
        last_visc = lederer_roegiers(self.params[0], normalise(moles[:2]), visc_parts[2:])
        last_mole = sum(moles[:2])
        partitions = self.components_amount - 2
        for i in range(partitions):
            last_visc = lederer_roegiers(self.params[i + 1])


def main_visc_check():
    substance = "R115"
    data = ViscData(substance)
    pref = 0.06

    visc = data.visc_d
    t = [258, 273.16, 288.14, 303.13, 318.11, 333.1, 348.08, 363.26]
    t = [i - 273.15 for i in t]
    # oil_visc_d = [412.33, 136.68, 57.08, 28.17, 16.01, 9.83, 6.69, 4.73]
    # visc = {0.06: {t[i]: oil_visc_d[i] for i in range(len(t))}}
    temps = list(visc[pref].keys())
    viscs = [i[2] for i in data.get_by_pressure(2, pref)]

    test = VFT()
    test.setup(pref, visc)
    print(
        f'D = {test.params[0]}\nE = {test.params[1]}\nF = {test.params[2]}\nG = {test.params[3:6]}\nH = {test.params[6]}')

    font = {'family': 'Times New Roman', 'size': 10.5}
    matplotlib.rc('font', **font)
    fig, ax = plt.subplots()
    ax.scatter(temps, viscs, label="Экспериментальные данные")
    viscs_calc = [test.calc_visc(i, pref) for i in arange(temps[0], temps[-1] + 1)]

    viscs_temp = [test.calc_visc_def(i) for i in arange(temps[0], temps[-1] + 1)]
    ax.plot(arange(temps[0], temps[-1] + 1), viscs_calc, "r",
            label="Расчётные значения")
    ax.set_ylabel("$\it{Вязкость, сантипуаз}$")
    ax.set_xlabel('$\it{T, °C}$')
    plt.grid()
    plt.legend()
    plt.show()


def mix_visc_check():
    data_r115 = ViscData("R115")
    data_r22 = ViscData("R22")
    pref = 0.06
    t_list = [-10, 0, 20, 40, 70, 80, 125]
    display_t = 125
    step = 0.1
    '''visc_r115 = data_r115.visc_d
    visc_r22 = data_r22.visc_d

    visc_oil = {pref: {-15.171: 413.005, -10.1111: 277.596, -5.0512: 191.427, 0.00868263: 137.18, 4.87395: 99.5737,
                       9.93383: 75.1096, 14.9937: 57.3869, 20.0536: 44.4115, 24.9189: 35.2622, 29.9787: 28.359,
                       35.0386: 23.1013, 39.9039: 19.0611, 44.9638: 15.9304, 50.0237: 13.4856, 54.8889: 11.5632,
                       59.9488: 9.91492, 65.0087: 8.7223, 70.0686: 7.67312, 74.9338: 6.75015, 79.9937: 6.01479,
                       85.0536: 5.35955, 90.1135: 4.83728, 95.1734: 4.31031, 100.039: 3.94046}}'''
    # Viscosity solvers for oil, r115 and r22
    vft_oil = VFT()
    vft_oil.params = [0.0543, 940.01, 152.78, 0, 0, 0, 0]
    vft_r115 = VFT()
    vft_r115.params = [9.0340e-05, -1848.5, 1174.4, 0, 0, 0, 0]
    vft_r22 = VFT()
    vft_r22.params = [0.0186, 674.46, 391.24, 0, 0, 0, 0]


    density_oil = 
    mole_arr = []

    '''print('Initializing oil setup')
    if vft_oil.setup_def():
        print('vft_oil setup complete')
    else:
        raise Exception('vft_oil setup incomplete, aborting')
    print('Initializing vft_r115 setup')
    if vft_r115.setup():
        print('vft_r115 setup complete')
    else:
        raise Exception('vft_r115 setup incomplete, aborting')
    print('Initializing vft_r22 setup')
    if vft_r22.setup():
        print('vft_r22 setup complete')
    else:
        raise Exception('vft_r22 setup incomplete, aborting')'''
    print('Initializing mole frac solver')
    mole_frac_solver = create_frac_calc()
    print('mole frac solver setup complete')

    print('Loading cavestri data for mixture')
    cavestri_data = get_data(True)
    print('Data loading complete')
    points_mu_t = cavestri_data.get_by_t(display_t)
    p_scatter = [i.p for i in points_mu_t]
    mu_scatter = [i.mu for i in points_mu_t]

    print('Initializing calculation of viscosity points')
    p_calc = arange(min(p_scatter), max(p_scatter) + step, step)
    visc_calc = []
    for p in p_calc:
        mole_frac = mole_frac_solver.calculate_flash(display_t, p)[0]
        mu_r22_p = vft_r22.calc_visc_def(display_t)
        mu_r115_p = vft_r115.calc_visc_def(display_t)
        mu_oil_p = vft_oil.calc_visc_def(display_t)
        mu_list = [mu_r22_p, mu_r115_p]
        visc_refr = kendall_monroe(visc_list, normalise(mole_frac[:2]))
        visc_calc.append(kendall_monroe([visc_refr, visc_oil_p], [sum(mole_frac[:2]), mole_frac[2]]))
    print("Calculation complete, initiating plotting")
    # Displayed pressure
    p_plot = 0.06
    for t_cur in t_list:
        t_cur_points = cavestri_data.get_by_t(t_cur)
        points = [i.p for i in t_cur_points]
        x_exp = [i.mu for i in t_cur_points]
        x_calc = []
        for p in points:
            mole_frac = mole_frac_solver.calculate_flash(t_cur, p)[0]
            visc_r22_p = vft_r22.calc_visc_def(t_cur)
            visc_r115_p = vft_r115.calc_visc_def(t_cur)
            visc_oil_p = vft_oil.calc_visc_def(t_cur)
            visc_list = [visc_r22_p, visc_r115_p]
            visc_refr = kendall_monroe(visc_list, normalise(mole_frac[:2]))
            x_calc.append(kendall_monroe([visc_refr, visc_oil_p], [sum(mole_frac[:2]), mole_frac[2]]))
        print(f'AAD для t={t_cur}:\t{aad(x_exp, x_calc)}')

    font = {'family': 'Times New Roman', 'size': 10.5}
    matplotlib.rc('font', **font)
    fig, ax = plt.subplots()
    ax.scatter(p_scatter, mu_scatter, label="Экспериментальные данные")
    ax.plot(p_calc, visc_calc, "r", label="Расчётные значения")
    # ax.plot(arange(temps[0], temps[-1] + 1), viscs_calc, "r",
    #        label="Расчётные значения")
    ax.set_ylabel("$\it{Вязкость, сантипуаз}$")
    ax.set_xlabel('$\it{P, МПа}$')
    plt.grid()
    plt.legend()
    plt.show()


if __name__ == "__main__":
    mix_visc_check()
