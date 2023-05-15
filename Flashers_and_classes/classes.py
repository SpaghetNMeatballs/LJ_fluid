from functions.read_func import read_data_from_file_lines
from thermo import ChemicalConstantsPackage
from scipy.optimize import minimize
from thermo import CEOSGas, CEOSLiquid, FlashVL, GCEOSMIX
from dataclasses import dataclass


class FlashCalculator:
    def __init__(self, class_flash, properties, constants, initial_guess):
        self.class_flash = class_flash
        self.properties = properties
        self.constants = constants
        self.params = initial_guess

    def func_opt(self, params: list[float], fraction_list: list[list[list[float]]], t: list[float],
                 p: list[list[float]], zs: list[float]) -> float:
        count = 0
        sum1 = 0
        sum2 = 0
        fractions1 = fraction_list[0]
        fractions2 = fraction_list[1]
        for i in range(len(t)):
            for j in range(len(p[i])):
                k13, k23, k12 = [params[ind] + (t[i] - t[0]) * params[ind + 1] for ind in
                                 range(0, 5, 2)]
                kijs = [[0, k12, k13],
                        [k12, 0, k23],
                        [k13, k23, 0]]

                eos_kwargs = {'Pcs': self.constants.Pcs, 'Tcs': self.constants.Tcs, 'omegas': self.constants.omegas,
                              'kijs': kijs}
                gas = CEOSGas(self.class_flash, eos_kwargs=eos_kwargs)
                liquid = CEOSLiquid(self.class_flash, eos_kwargs=eos_kwargs)
                flasher = FlashVL(self.constants, self.properties, liquid=liquid, gas=gas)

                flashed = flasher.flash(T=t[i] + 273.15, P=p[i][j] * pow(10, 5), zs=zs)
                if flashed.VF != 1:
                    sum1 += abs(fractions1[i][j] - flashed.liquid0.zs[0])
                    sum2 += abs(fractions2[i][j] - flashed.liquid0.zs[1])
                else:
                    sum1 += abs(fractions1[i][j] - 0.0)
                    sum2 += abs(fractions2[i][j] - 0.0)
                count += 1

        avg_error1 = sum1 / count
        avg_error2 = sum2 / count
        return avg_error1 + avg_error2

    def optimize_parameters(self, initial_guess_k13_k23: list[float], fractions: callable, t: list[float],
                            p: list[list[float]], zs: list[float]):
        args = (fractions, t, p, zs)
        bounds = ((-1, 1), (0, 0.004), (-1, 1), (0, 0.004), (-1, 1), (0, 0.004))
        result = minimize(self.func_opt, initial_guess_k13_k23, args=args, bounds=bounds)
        if result.success:
            self.params = result.x
        else:
            raise ValueError(result.message)
        return result

    def calculate_flash(self, t: list[float], p: list[list[float]], fractions: list[float], t_ind: int, p_ind: int):
        k13 = self.params[0] + (t[t_ind] - t[0]) * self.params[1]
        k23 = self.params[2] + (t[t_ind] - t[0]) * self.params[3]
        k12 = self.params[4] + (t[t_ind] - t[0]) * self.params[5]

        kijs = [[0, k12, k13],
                [k12, 0, k23],
                [k13, k23, 0]]

        eos_kwargs = {'Pcs': self.constants.Pcs, 'Tcs': self.constants.Tcs, 'omegas': self.constants.omegas,
                      'kijs': kijs}
        gas = CEOSGas(self.class_flash, eos_kwargs=eos_kwargs)
        liquid = CEOSLiquid(self.class_flash, eos_kwargs=eos_kwargs)
        flasher = FlashVL(self.constants, self.properties, liquid=liquid, gas=gas)

        PT = flasher.flash(T=t[t_ind] + 273.15, P=p[t_ind][p_ind] * pow(10, 5), zs=fractions)

        if PT.VF == 1:  # Есть только пар
            return [PT.gas.zs, t, p, PT.VF, k12, k13, k23]
        else:  # Не только пар
            return [PT.liquid0.zs, t, p, PT.VF, k12, k13, k23]


@dataclass
class Oil:
    M_oil: float
    T_oil: float
    P_oil: float
    Omega_oil: float


def create_constants_package(oil: Oil, IDs) -> ChemicalConstantsPackage:
    # Перезаписали критические параметры для воды под масло ISO 32
    constants = ChemicalConstantsPackage.from_IDs(IDs)[0]
    return ChemicalConstantsPackage(
        MWs=constants.MWs[:-1] + [oil.M_oil],
        names=constants.names[:-1] + ['ISO 32'],
        omegas=constants.omegas[:-1] + [oil.Omega_oil],
        Pcs=constants.Pcs[:-1] + [oil.P_oil],
        Tcs=constants.Tcs[:-1] + [oil.T_oil]
    )


'''

# Плотность R-22:
oil_PRMIX = Oil(400.0, 800, 950e3, 0.8)

# Свойства Oil (Молярная масса масла, температура (Кельвины), давление (Па), константа омега)
oil_SRKMIX = Oil(400.0, 800, 950e3, 0.2)

# Создание констант и экземпляров классов
IDs = ['r-22', 'chloropentafluoroethane', 'water']

chemical_constants_PR = create_constants_package(oil_PRMIX, IDs)

chemical_constants_SR = create_constants_package(oil_SRKMIX, IDs)

temperature = read_data_from_file_lines('temperature.txt')
pressure = read_data_from_file_lines('pressure.txt')
viscosity = read_data_from_file_lines('viscosity.txt')
density = read_data_from_file_lines('density.txt')
Fraction_R22 = read_data_from_file_lines('fraction_r22.txt')
Fraction_R115 = read_data_from_file_lines('fraction_r115.txt')
'''
