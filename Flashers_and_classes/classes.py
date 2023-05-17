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
        if len(fraction_list) == 2:
            fractions2 = fraction_list[1]
        for i in range(len(t)):
            for j in range(len(p[i])):
                if len(zs) == 3:
                    k13, k23, k12 = [params[ind] + (t[i] - t[0]) * params[ind + 1] for ind in
                                     range(0, 5, 2)]
                    kijs = [[0, k12, k13],
                            [k12, 0, k23],
                            [k13, k23, 0]]
                elif len(zs) == 2:
                    k12 = params[0] + (t[i] - t[0]) * params[1]
                    kijs = [[0, k12],
                            [k12, 0]]

                eos_kwargs = {'Pcs': self.constants.Pcs, 'Tcs': self.constants.Tcs, 'omegas': self.constants.omegas,
                              'kijs': kijs}
                gas = CEOSGas(self.class_flash, eos_kwargs=eos_kwargs)
                liquid = CEOSLiquid(self.class_flash, eos_kwargs=eos_kwargs)
                flasher = FlashVL(self.constants, self.properties, liquid=liquid, gas=gas)

                flashed = flasher.flash(T=t[i] + 273.15, P=p[i][j] * pow(10, 5), zs=zs)
                if flashed.VF != 1:
                    sum1 += abs(fractions1[i][j] - flashed.liquid0.zs[0])
                    if len(zs) == 3:
                        sum2 += abs(fractions2[i][j] - flashed.liquid0.zs[1])
                else:
                    sum1 += abs(fractions1[i][j] - 0.0)
                    if len(zs) == 3:
                        sum2 += abs(fractions2[i][j] - 0.0)
                count += 1

        avg_error1 = sum1 / count
        if len(zs) == 3:
            avg_error2 = sum2 / count
        if len(zs) == 3:
            return avg_error1 + avg_error2
        elif len(zs) == 2:
            return avg_error1

    def optimize_parameters(self, initial_guess_k13_k23: list[float], fractions: callable, t: list[float],
                            p: list[list[float]], zs: list[float]):
        args = (fractions, t, p, zs)
        bounds_zs_3 = ((-1, 1), (0, 0.004), (-1, 1), (0, 0.004), (-1, 1), (0, 0.004))
        bounds_zs_2 = ((-1, 1), (0, 0.004))
        if len(zs) == 3:
            result = minimize(self.func_opt, initial_guess_k13_k23, args=args, bounds=bounds_zs_3)
        elif len(zs) == 2:
            result = minimize(self.func_opt, initial_guess_k13_k23, args=args, bounds=bounds_zs_2)
        if result.success:
            self.params = result.x
        else:
            raise ValueError(result.message)
        return result

    def calculate_flash(self, t: list[float], p: list[list[float]], fractions: list[float], t_ind: int, p_ind: int):
        if len(fractions) == 3:
            k13 = self.params[0] + (t[t_ind] - t[0]) * self.params[1]
            k23 = self.params[2] + (t[t_ind] - t[0]) * self.params[3]
            k12 = self.params[4] + (t[t_ind] - t[0]) * self.params[5]

            kijs = [[0, k12, k13],
                    [k12, 0, k23],
                    [k13, k23, 0]]
        elif len(fractions) == 2:
            k12 = self.params[0] + (t[t_ind] - t[0]) * self.params[1]
            kijs = [[0, k12],
                    [k12, 0]]

        eos_kwargs = {'Pcs': self.constants.Pcs, 'Tcs': self.constants.Tcs, 'omegas': self.constants.omegas,
                      'kijs': kijs}
        gas = CEOSGas(self.class_flash, eos_kwargs=eos_kwargs)
        liquid = CEOSLiquid(self.class_flash, eos_kwargs=eos_kwargs)
        flasher = FlashVL(self.constants, self.properties, liquid=liquid, gas=gas)

        PT = flasher.flash(T=t[t_ind] + 273.15, P=p[t_ind][p_ind] * pow(10, 5), zs=fractions)

        if PT.VF == 1:  # Есть только пар
            return [PT.gas.zs, t[t_ind], p[t_ind][p_ind], PT.VF, k12] + ([k13, k23] if len(fractions) == 3 else [])
        else:  # Не только пар
            return [PT.liquid0.zs, t[t_ind], p[t_ind][p_ind], PT.VF, k12] + ([k13, k23] if len(fractions) == 3 else [])


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


