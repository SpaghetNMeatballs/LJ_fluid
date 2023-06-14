import matplotlib.pyplot as plt
from numpy import arange, array, zeros
from scipy.optimize import minimize
import matplotlib
from calc_functions import DEF, HG, grunberg_nissan, lederer_roegiers, arrhenius, kendall_monroe, refutas
from common import Point, ViscData, aad, aad_new, create_frac_calc, create_frac_calc_2comp, normalise, mole_to_mass
from data_csv.cavestri_loader import get_data
from Flashers_and_classes.classes import FlashCalculator
from density import Density


class VFT:
    def __init__(self):
        self.pref = 0
        self.params = [0., 0., 0., 0., 0., 0., 0.]
        self.maxiter = 2 ** 20

    def setup_def(self, mu_pref: dict[float:float]) -> bool:
        visc_pref = list(mu_pref.values())
        temps = list(mu_pref.keys())

        # Выбираем начальное приближение на основе характера экспериментальных данных
        x0 = [0.000852102519255135, -67.82280822979249, 4.052812591576497]

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
        def_flag = self.setup_def(mu_data[pref])
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
        f'D = {test.params[0]}\nE = {test.params[1]}\nF = {test.params[2]}\nG = {test.params[3:6]}\nH = {test.params[6]}')

    font = {'family': 'Times New Roman', 'size': 14}
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


def oil_coeff_calc():
    t_arr = [-25, 0, 25, 50, 75, 100, 125]
    t_calc = arange(t_arr[0], t_arr[-1] + 0.1, 0.1)
    rho_arr = [1733.38, 289.693, 64.8321, 19.9096, 7.79992, 4.19357, 3.09385]
    oil_calc = VFT()
    oil_calc.setup_def({t_arr[i]: rho_arr[i] for i in range(len(t_arr))})
    print(oil_calc.params)
    font = {'family': 'Times New Roman', 'size': 14}
    matplotlib.rc('font', **font)
    fig, ax = plt.subplots()
    ax.scatter(t_arr, rho_arr, label="Экспериментальные данные")
    ax.plot(t_calc, [oil_calc.calc_visc_def(i) for i in t_calc], "r", label="Расчётные значения")
    ax.set_ylabel("$\it{Вязкость, сантипуаз}$")
    ax.set_xlabel('$\it{P, бар}$')
    plt.yscale('log')
    plt.grid()
    plt.legend()
    plt.show()


def mix_visc_check():
    t_list = [-10, 0, 20, 40, 70, 80, 125]
    display_t = 20
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
    vft_oil.params = [0.005337316178523319, 1755.4979644152024, 109.82206931834386, 0, 0, 0, 0]
    vft_r115 = VFT()
    vft_r115.params = [0.0008577907640623858, -69.08004012663991, 4.2625169264404175, 0, 0, 0, 0]
    vft_r22 = VFT()
    vft_r22.params = [2.2755471578933423e-05, 626.8283963685842, 2.50475798005479, 0, 0, 0, 0]

    point0_r22 = Point(t=-60, p=0.07, rho=1463.67431328749)
    point0_r115 = Point(t=-49.5, p=0.07, rho=1581.88691446926)
    point0_oil = Point(t=-25, rho=0.935)

    density_r22 = Density(point0_r22)
    density_r115 = Density(point0_r115)
    density_oil = Density(point0_oil)

    # Присваиваем параматреы, расчитанные ранее
    density_r22.e, density_r22.beta = 132772382183158.875000, 0.001752
    density_r115.e, density_r115.beta = 21106748519192.570312, 0.001916
    density_oil.beta = 0.0010332054833384788

    print('Initializing mole frac solver')
    mole_frac_solver = create_frac_calc_2comp()
    print('Mole frac solver setup complete')

    print('Loading cavestri data for mixture')
    cavestri_data = get_data(False)
    print('Data loading complete')
    points_mu_t = cavestri_data.get_by_t(display_t)
    p_scatter = [i.p for i in points_mu_t]
    mu_scatter = [i.mu for i in points_mu_t]

    print('Initializing calculation of viscosity points')

    p_calc = arange(min(p_scatter), max(p_scatter) + step, step)
    visc_calc = []
    points_minimize = []
    exp_minimize = []
    for i in cavestri_data.points:
        mole_frac = mole_frac_solver.calculate_flash(i.t, i.p)[0]

        mu_r22 = vft_r22.calc_visc_def(i.t)

        v_oil_p = vft_oil.calc_visc_def(i.t)
        rho_oil_p = density_oil.calc_rho_t(i.t)
        mu_oil = v_oil_p * rho_oil_p

        points_minimize.append([mu_r22, mu_oil, mole_frac[0], mole_frac[1]])
        exp_minimize.append(i.mu)

    args_lr = (lederer_roegiers, points_minimize, exp_minimize)
    args_gn = (grunberg_nissan, points_minimize, exp_minimize)
    gn_minimize_result = minimize(aad_new,
                                  x0=[9.468922609917012, 7.457512076395755, 5.218073731487143, 1.2839020316809933],
                                  method='Nelder-Mead', args=args_gn, tol=1e-8,
                                  options={'maxiter': 2 ** 20, 'maxfev': 2 ** 20})
    lr_minimize_result = minimize(aad_new, x0=10.430192947387695, method='Nelder-Mead', args=args_lr, tol=1e-5,
                                  options={'maxiter': 2 ** 20, 'maxfev': 2 ** 20})
    alpha = lr_minimize_result.x
    params_gn = gn_minimize_result.x

    for p in p_calc:
        cur_point = Point(t=display_t, p=p * 10)
        mole_frac = mole_frac_solver.calculate_flash(display_t, p)[0]
        mass_frac = mole_to_mass(mole_frac, [86.47, 400.0])

        mu_r22_p = vft_r22.calc_visc_def(display_t)
        rho_r22_p = density_r22.calc_rho_full(cur_point)
        v_r22_p = mu_r22_p / rho_r22_p * 10 ** 6

        '''mu_r115_p = vft_r115.calc_visc_def(display_t)
        rho_r115_p = density_r115.calc_rho_full(cur_point)
        v_r115_p = mu_r115_p / rho_r115_p'''

        v_oil_p = vft_oil.calc_visc_def(display_t)
        rho_oil_p = density_oil.calc_rho_t(display_t)
        mu_oil_p = v_oil_p * rho_oil_p

        # mu_refr_list = [mu_r22_p, mu_r115_p]
        # mu_refr = arrhenius(lr_params[0], mu_refr_list, normalise(mole_frac[:2]))
        visc_calc.append(grunberg_nissan(params_gn, [mu_r22_p, mu_oil_p, mole_frac[0], mole_frac[1]]))
        # visc_calc.append(lederer_roegiers(alpha, [mu_r22_p, mu_oil_p, mole_frac[0], mole_frac[1]]))
    print("Calculation complete, initiating plotting")

    for t_cur in t_list:
        t_cur_points = cavestri_data.get_by_t(t_cur)
        if len(t_cur_points) == 0:
            continue
        points = [i.p for i in t_cur_points]
        x_exp = [i.mu for i in t_cur_points]
        x_calc = []
        for i in t_cur_points:
            mole_frac = mole_frac_solver.calculate_flash(t_cur, i.p)[0]
            mass_frac = mole_to_mass(mole_frac, [86.47, 400.0])

            mu_r22_p = vft_r22.calc_visc_def(t_cur)
            rho_r22_p = density_r22.calc_rho_full(i)
            v_r22_p = mu_r22_p / rho_r22_p

            # mu_r115_p = vft_r115.calc_visc_def(display_t)
            # rho_r115_p = density_r115.calc_rho_full(cur_point)
            # v_r115_p = mu_r115_p / rho_r115_p

            v_oil_p = vft_oil.calc_visc_def(t_cur)
            rho_oil_p = density_oil.calc_rho_t(t_cur)
            mu_oil_p = v_oil_p * rho_oil_p

            # mu_refr_list = [mu_r22_p, mu_r115_p]
            # mu_refr = arrhenius(lr_params[0], mu_refr_list, normalise(mole_frac[:2]))
            # visc_calc.append(arrhenius(lr_params[1], [mu_refr, mu_oil_p], [sum(mole_frac[:2]), mole_frac[2]]))
            x_calc.append(grunberg_nissan(params_gn, [mu_r22_p, mu_oil_p, mole_frac[0], mole_frac[1]]))
        print(f'AAD для t={t_cur}:\t{aad(x_exp, x_calc)}')

    font = {'family': 'Times New Roman', 'size': 14}
    matplotlib.rc('font', **font)
    fig, ax = plt.subplots()
    ax.scatter(p_scatter, mu_scatter, label="Экспериментальные данные")
    ax.plot(p_calc, visc_calc, "r", label="Расчётные значения")
    # ax.plot(arange(temps[0], temps[-1] + 1), viscs_calc, "r",
    #        label="Расчётные значения")
    ax.set_ylabel("$\it{Вязкость, сантипуаз}$")
    ax.set_xlabel('$\it{P, бар}$')
    plt.grid()
    plt.legend()
    plt.show()


if __name__ == "__main__":
    mix_visc_check()
