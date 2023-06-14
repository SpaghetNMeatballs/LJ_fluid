from common import Point, aad_new, ViscData
from calc_functions import rho_t, rho_p, rho_full, rho_paredes_t, rho_paredes
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib
from numpy import arange


class Density:
    def __init__(self, point0: Point):
        self.point0 = point0
        self.beta = 0
        self.e = 6.089e13

    def calc_rho_t(self, t: float) -> float:
        return rho_t(self.beta, [self.point0.rho, self.point0.t], t)

    def calc_rho_p(self, p: float) -> float:
        return rho_p(self.beta, [self.point0.rho, self.point0.p], p)

    def calc_rho_full(self, point: Point) -> float:
        p = point.p
        t = point.t
        return rho_full([self.e, self.beta], [self.point0.rho, self.point0.p, self.point0.t], [p, t])

    def setup_rho_t(self, points: list[Point], p0_ind: int) -> bool:
        if self.point0:
            assert self.point0 == points[p0_ind], "Wrong data passed"
        self.point0 = points[p0_ind]
        temps = [i.t for i in points]
        rho = [i.rho for i in points]
        args = (rho_t, temps, rho, [self.point0.rho, self.point0.t])
        result = minimize(fun=aad_new, x0=0, method='Powell', args=args)
        self.beta = result.x[0]
        print(result.message)
        return result.success

    def setup_rho_p(self, p: list[float], rho: list[float]) -> bool:
        args = (rho_p, p, rho, [self.point0.rho, self.point0.p])
        result = minimize(fun=aad_new, x0=3, args=args)
        self.e = result.x[0]
        return result.success

    def setup_rho_full(self, rho, points) -> bool:
        static_params = (self.point0.rho, self.point0.p, self.point0.t)
        args = (rho_full, points, rho, static_params)
        result = minimize(fun=aad_new, method='Nelder-Mead', x0=[3, 0], args=args,
                          options={'maxiter': 100000, 'maxfev': 100000})
        self.e, self.beta = result.x
        print(result.message)
        return result.success


class DensityParedes:
    def __init__(self, pref: float):
        self.pref = pref
        self.params = [1232.7, -0.8404, 0.0001143, 0.08365, 452.16, -1.5087, 0.0014034]

    def setup_t(self, rho_pref: dict[float:float]) -> bool:
        t = list(rho_pref.keys())
        rho = list(rho_pref.values())
        args = (rho_paredes_t, t, rho)
        minimized = minimize(fun=aad_new, method='Nelder-Mead', x0=[354, -7, 0], args=args,
                             options={'maxiter': 100000, 'maxfev': 100000})
        print(minimized.message)
        print(f'AAD_t = {aad_new(minimized.x, rho_paredes_t, t, rho)}')
        self.params[:3] = minimized.x
        return minimized.success

    def setup(self, rho_dict: dict[float:dict[float:float]], pref: float = None) -> bool:
        p_list = list(rho_dict.keys())
        if pref is None:
            pref = p_list[int(len(p_list) / 2)]
        if not self.setup_t(rho_dict[pref]):
            raise Exception("Error while optimizing t equation")
        points = []
        rho = []
        for i in p_list:
            t_list = list(rho_dict[i].keys())
            for j in t_list:
                points.append([i, j])
                rho.append(rho_dict[i][j])
        args = (rho_paredes, points, rho, self.params[:3] + [pref])
        x0 = [0.08365, 452.16, -1.5087, 0.0014034]
        minimized = minimize(fun=aad_new, x0=x0, args=args,
                             options={'maxiter': 100000, 'maxfev': 100000})
        print(minimized.message)
        print(f'AAD = {aad_new(minimized.x, rho_paredes, points, rho, self.params[:3] + [pref])}')
        self.params[3:] = minimized.x
        return minimized.success

    def calc_t(self, t: float) -> float:
        return rho_paredes_t(self.params[:3], t)


def rho_test():
    substance = "R22"
    substance_data = ViscData(substance)
    points_to_process = substance_data.points
    points = [[i.p, i.t] for i in points_to_process]
    rho = [i.rho for i in points_to_process]
    point0 = points_to_process[int(len(points_to_process) / 2)]
    density_calc = Density(point0)
    if not density_calc.setup_rho_full(rho, points):
        print('Error while optimizing')

    print(f'rho0 = {density_calc.point0.rho}\nt0 = {density_calc.point0.t}\np0 = {density_calc.point0.p}')
    print("Beta = %f" % density_calc.beta)
    print("E = %f" % density_calc.e)
    static_params = (density_calc.point0.rho, density_calc.point0.p, density_calc.point0.t)
    t_arr = [i.t for i in substance_data.get_by_pressure(point0.p)]
    rho_arr = [i.rho for i in substance_data.get_by_pressure(point0.p)]
    t_calc = arange(min(t_arr), max(t_arr) + 0.1, 0.1)
    print(f'AAD = {aad_new([3, 0], rho_full, points, rho, static_params)}%')
    font = {'family': 'Times New Roman', 'size': 14}
    matplotlib.rc('font', **font)
    fig, ax = plt.subplots()
    ax.scatter(t_arr, rho_arr, label="Экспериментальные данные")
    ax.plot(t_calc, [density_calc.calc_rho_full(Point(t=i, p=point0.p)) for i in t_calc], "r",
            label="Расчётные значения")
    ax.set_ylabel("Rho, м3/кг")
    ax.set_xlabel("P, МПа")
    plt.grid()
    plt.legend()
    plt.show()


def rho_oil_tait():
    t_arr = [-25, 0, 25, 50, 75, 100, 125]
    rho_arr = [1733.38, 289.693, 64.8321, 19.9096, 7.79992, 4.19357, 3.09385]
    t_calc = arange(t_arr[0], t_arr[-1] + 0.1, 0.1)
    rho_calc = DensityParedes(0)
    rho_calc.setup_t({t_arr[i]: rho_arr[i] for i in range(len(t_arr))})
    fig, ax = plt.subplots()
    ax.scatter(t_arr, rho_arr, label="Экспериментальные данные")
    ax.plot(t_calc, [rho_calc.calc_t(i) for i in t_calc], "r", label="Расчётные значения")
    ax.set_ylabel("Rho, м3/кг")
    ax.set_xlabel("P, МПа")
    plt.legend()
    plt.show()


def rho_oil():
    t_arr = [-25, 0, 25, 50, 75, 100, 125]
    rho_arr = [i * 10 ** 3 for i in [0.935, 0.910277, 0.888055, 0.86833, 0.847776, 0.828056, 0.807502]]
    point_list = [Point(t=t_arr[i], rho=rho_arr[i]) for i in range(len(t_arr))]
    t_calc = arange(t_arr[0], t_arr[-1] + 0.1, 0.1)
    rho_calc = Density(point_list[4])
    rho_calc.setup_rho_t(point_list, 4)
    print(rho_calc.beta)
    print(aad_new(rho_calc.beta, rho_t, t_arr, rho_arr, [rho_calc.point0.rho, rho_calc.point0.t]))
    font = {'family': 'Times New Roman', 'size': 14}
    matplotlib.rc('font', **font)
    fig, ax = plt.subplots()
    ax.scatter(t_arr, rho_arr, label="Экспериментальные данные")
    ax.plot(t_calc, [rho_calc.calc_rho_t(i) for i in t_calc], "r", label="Расчётные значения")
    ax.set_ylabel("Rho, м3/кг")
    plt.grid()
    ax.set_xlabel("P, МПа")
    plt.legend()
    plt.show()


def rho_paredes_test():
    substance = "R22"
    substance_data = ViscData(substance)
    density_dict = {}
    for i in substance_data.points:
        if i.p in list(density_dict.keys()):
            density_dict[i.p][i.t] = i.rho
        else:
            density_dict[i.p] = {i.t: i.rho}
    pref = 0.06
    rho_pref = density_dict[pref]

    test = DensityParedes(pref)
    if not test.setup(density_dict, pref):
        raise Exception('Error while optimizing')

    pcalc = 0.06
    fig, ax = plt.subplots(1, 1)
    t_exp = list(density_dict[pcalc].keys())
    rho_exp = list(density_dict[pcalc].values())
    t_calc = arange(t_exp[0], t_exp[-1] + 0.1, 0.1)

    font = {'family': 'Times New Roman', 'size': 10.5}
    matplotlib.rc('font', **font)
    ax.plot(t_calc, [test.calc_t(i) for i in t_calc], "r", label="Расчётные значения")
    ax.scatter(t_exp, rho_exp, label="Экспериментальные данные")
    ax.set_ylabel("$\it{Плотность, кг/м3}$")
    ax.set_xlabel('$\it{T, °C}$')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    rho_oil()
