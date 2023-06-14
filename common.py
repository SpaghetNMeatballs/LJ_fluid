import os
from dataclasses import dataclass
from typing import Optional
from Flashers_and_classes.classes import FlashCalculator, create_constants_package, Oil
from functions.read_func import read_data_from_file_lines
from functions.chemical_functions import mass_fractions_to_mole_fractions, calculate_acentrick_factor
from thermo import SRKMIX
from thermo.chemical import Chemical
from math import log


@dataclass
class Point:
    t: Optional[float] = None
    p: Optional[float] = None
    x: Optional[list[float]] = None
    mu: Optional[float] = None
    v: Optional[float] = None
    rho: Optional[float] = None

    def __eq__(self, other):
        if self.t != other.t:
            return False
        if self.p != other.p:
            return False
        if self.mu != other.mu:
            return False
        if self.rho != other.rho:
            return False
        x_size = -1
        o_size = -1
        if self.x:
            x_size = len(self.x)
        if other.x:
            o_size = len(other.x)
        if x_size != o_size:
            return False
        for i in range(x_size):
            if self.x[i] != other.x[i]:
                return False
        return True


def mole_to_mass(mole_frac: list[float], mole_list: list[float]) -> list[float]:
    result = mole_frac.copy()
    result = [result[i] * mole_list[i] for i in range(len(result))]
    temp_sum = sum(result)
    result = [result[i] / temp_sum for i in range(len(result))]
    return result


def create_frac_calc_2comp() -> FlashCalculator:
    g = 0.55
    zs = [g, 1 - g]  # R-22, ISO 32
    T = read_data_from_file_lines('temperature_r22.txt')
    t_min = T[0]  # начальная температура для оптимизации
    P = read_data_from_file_lines('pressure_r22.txt')
    Fraction_R22 = read_data_from_file_lines('fraction_r22_alone.txt')
    IDs = ['r-22', 'water']

    # Молярная масса хладагента R-22
    r22 = Chemical(IDs[0])
    M_r22 = round(r22.MW, 2)  # 86.47

    # Критические параметры масла:
    M_oil = 400.0  # Молярная масса масла
    T_oil = 800  # Кельвины
    P_oil = 950e3  # Па
    # Расчет ацентрического фактора по корреляции Ли и Кесслера
    omega = calculate_acentrick_factor(T_oil, P_oil, M_oil)

    # Расчет мольных долей для каждого экспериментального компонента системы
    mole_fractions_r22, mole_fractions_oil = mass_fractions_to_mole_fractions(
        Fraction_R22, None, M_r22, None, M_oil
    )

    # Расчет мольных долей для каждого компонента системы реализуемой программы
    zs_r22_mole, zs_oil_mole = mass_fractions_to_mole_fractions(
        [[zs[0]]], None, M_r22, None, M_oil
    )
    zs_mole = [zs_r22_mole[0][0], zs_oil_mole[0][0]]

    oil_SRKMIX = Oil(M_oil, T_oil, P_oil, omega)
    chemical_constants_SR = create_constants_package(oil_SRKMIX, IDs)

    # ================= Оптимизация коэффициентов =================
    # Оптимизируемые коэффициенты для R-22
    initial_guess_k13_k23 = [-0.1, 0.001]
    flash_calculator = FlashCalculator(SRKMIX, oil_SRKMIX, chemical_constants_SR, initial_guess_k13_k23, zs_mole, t_min)
    result = flash_calculator.optimize_parameters(initial_guess_k13_k23, [mole_fractions_r22],
                                                  T, P, zs_mole)
    return flash_calculator


def create_frac_calc() -> FlashCalculator:
    g = 0.7
    zs = [0.488 * g, 0.512 * g, 1 - g]  # R-22, R-115, ISO 32
    T = read_data_from_file_lines('temperature.txt')
    t_min = T[0]  # начальная температура для оптимизации
    P = read_data_from_file_lines('pressure.txt')
    fractions = [read_data_from_file_lines('fraction_r22.txt'), read_data_from_file_lines('fraction_r115.txt')]
    Fraction_R22, Fraction_R115 = fractions
    IDs = ['r-22', 'chloropentafluoroethane', 'water']

    # Молярная масса хладагента R-22
    r22 = Chemical(IDs[0])
    M_r22 = round(r22.MW, 2)  # 86.47

    # Молярная масса хладагента R-115
    r115 = Chemical(IDs[1])
    M_r115 = round(r115.MW, 2)  # 154.47

    # Критические параметры масла:
    M_oil = 400.0  # Молярная масса масла
    T_oil = 800  # Кельвины
    P_oil = 950e3  # Па
    # Расчет ацентрического фактора по корреляции Ли и Кесслера
    omega = calculate_acentrick_factor(T_oil, P_oil, M_oil)

    # Расчет мольных долей для каждого экспериментального компонента системы
    mole_fractions_r22, mole_fractions_r115, mole_fractions_oil = mass_fractions_to_mole_fractions(
        Fraction_R22, Fraction_R115, M_r22, M_r115, M_oil
    )

    # Расчет мольных долей для каждого компонента системы реализуемой программы
    zs_r22_mole, zs_r115_mole, zs_oil_mole = mass_fractions_to_mole_fractions(
        [[zs[0]]], [[zs[1]]], M_r22, M_r115, M_oil
    )
    zs_mole = [zs_r22_mole[0][0], zs_r115_mole[0][0], zs_oil_mole[0][0]]

    oil_SRKMIX = Oil(M_oil, T_oil, P_oil, omega)
    chemical_constants_SR = create_constants_package(oil_SRKMIX, IDs)

    # ================= Оптимизация коэффициентов =================
    # Оптимизируемые коэффициенты для R-22 и R-115
    initial_guess_k13_k23 = [0.21, 0.001, 0.35, 0.001, 0.1, 0.001]
    flash_calculator = FlashCalculator(SRKMIX, oil_SRKMIX, chemical_constants_SR, initial_guess_k13_k23, zs_mole, t_min)
    result = flash_calculator.optimize_parameters(initial_guess_k13_k23, [mole_fractions_r22, mole_fractions_r115],
                                                  T, P, zs_mole)
    if result.success:
        return flash_calculator
    else:
        raise Exception('Error while optimizing parameters of flash calculator')


def load_from_directory(dir_path: str, column: int) -> dict[float:dict[float:float]]:
    result = {}
    file_path = os.path.join(os.path.dirname(__file__), dir_path)
    for file in os.listdir(file_path):
        with open(file_path + f'\\{file}', "r") as f:
            text = f.read()
            # Считываем давление из имени
            p = float(file[:-4])
            result[p] = {}
            for line in text.split('\n'):
                cur_line = line.split(' ')
                # Забираем значение из искомой колонки
                t = float(cur_line[0])
                val = float(cur_line[column])
                result[p][t] = val
    return result


class ViscData:
    def __init__(self, dir_path: str):
        self.points = []
        density = load_from_directory(dir_path, 2)
        visc_k = load_from_directory(dir_path, 3)
        visc_d = load_from_directory(dir_path, 4)
        p = list(density.keys())
        for p_cur in p:
            for t_cur in list(density[p_cur].keys()):
                mu = visc_d[p_cur][t_cur]
                v = visc_k[p_cur][t_cur]
                rho = density[p_cur][t_cur]
                self.points.append(Point(t=t_cur, p=p_cur, v=v, mu=mu, rho=rho))

    def get_by_pressure(self, p: float) -> list[Point]:
        result = []
        for i in self.points:
            if i.p == p:
                result.append(i)
        return result

    def get_unique_pressure(self) -> list[float]:
        result = set()
        for i in self.points:
            result.add(i.p)
        result = list(result)
        result.sort()
        return result


def aad(x_exp, x_calc):
    assert (len(x_exp) == len(x_calc))
    returnable = 0
    for i in range(len(x_exp)):
        returnable += abs((x_exp[i] - x_calc[i]) / x_exp[i])
    returnable *= 100 / len(x_exp)
    return returnable


def aad_new(params: list[float], func: callable, points: list, exp: list[float],
            static_params: list[float] = None) -> float:
    if len(points) != len(exp):
        raise ValueError(f"Input lists have different shapes: {len(points)} and {len(exp)}")
    result = 0
    calc = []
    for i in range(len(points)):
        if static_params is None:
            calc.append(func(params, points[i]))

        else:
            calc.append(func(params, static_params, points[i]))

    for i in range(len(points)):
        result += abs((log(exp[i]) - log(calc[i])))

    result *= 100 / len(points)
    return result


def normalise(to_process: list[float]) -> list[float]:
    summ = sum(to_process)
    result = [i / summ for i in to_process]
    return result


if __name__ == "__main__":
    test_calc = create_frac_calc()
    print(test_calc.calculate_flash(-10, 3.38)[0])
