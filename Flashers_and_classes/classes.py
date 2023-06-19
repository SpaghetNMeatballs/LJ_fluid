from thermo import ChemicalConstantsPackage
from scipy.optimize import minimize
from thermo import CEOSGas, CEOSLiquid, FlashVL
from dataclasses import dataclass

'''
class FlashCalculator:
    def __init__(self, class_flash, properties, constants, initial_guess, fractions, t_min):
        self.class_flash = class_flash
        self.properties = properties
        self.constants = constants
        self.params = initial_guess
        self.fractions = fractions
        self.t_min = t_min

    def func_opt(self, params: list[float], fraction_list: list[list[list[float]]], t: list[float],
                 p: list[list[float]], zs: list[float]) -> float:
        count = 0
        sum1 = 0
        sum2 = 0
        fractions1 = fraction_list[0]
        if len(zs) == 3:
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

    def calculate_flash(self, t: float, p: float):
        if len(self.fractions) == 3:
            k13 = self.params[0] + (t - self.t_min) * self.params[1]
            k23 = self.params[2] + (t - self.t_min) * self.params[3]
            k12 = self.params[4] + (t - self.t_min) * self.params[5]

            kijs = [[0, k12, k13],
                    [k12, 0, k23],
                    [k13, k23, 0]]
        elif len(self.fractions) == 2:
            k12 = self.params[0] + (t - self.t_min) * self.params[1]
            kijs = [[0, k12],
                    [k12, 0]]

        eos_kwargs = {'Pcs': self.constants.Pcs, 'Tcs': self.constants.Tcs, 'omegas': self.constants.omegas,
                      'kijs': kijs}
        gas = CEOSGas(self.class_flash, eos_kwargs=eos_kwargs)
        liquid = CEOSLiquid(self.class_flash, eos_kwargs=eos_kwargs)
        flasher = FlashVL(self.constants, self.properties, liquid=liquid, gas=gas)

        PT = flasher.flash(T=t + 273.15, P=p * pow(10, 5), zs=self.fractions)

        if PT.VF == 1:  # Есть только пар
            return [PT.gas.zs, t, p, PT.VF, k12] + ([k13, k23] if len(self.fractions) == 3 else [])
        else:  # Не только пар
            return [PT.liquid0.zs, t, p, PT.VF, k12] + ([k13, k23] if len(self.fractions) == 3 else [])
'''


class FlashCalculator:
    def __init__(self, class_flash, properties, constants, initial_guess, fractions, t_min):
        self.class_flash = class_flash
        self.properties = properties
        self.constants = constants
        self.params = initial_guess
        self.fractions = fractions
        self.t_min = t_min

    def func_opt(self, params: list[float], fractions: dict[str, list[list[float]]], t: list[float],
                 p: list[list[float]], zs: list[float]) -> float:
        count = 0
        # Инициализируем листы в fractions
        gas_keys = list(fractions.keys())
        sum_dict = {key: 0 for key in gas_keys}

        for i in range(len(t)):
            for j in range(len(p[i])):
                kijs = [[0] * len(zs) for _ in range(len(zs))]
                index = 0
                for k in range(len(zs)):
                    for l in range(k + 1, len(zs)):
                        kijs[k][l] = kijs[l][k] = params[index] + (t[i] - t[0]) * params[index + 1]
                        index += 2

                eos_kwargs = {'Pcs': self.constants.Pcs, 'Tcs': self.constants.Tcs, 'omegas': self.constants.omegas,
                              'kijs': kijs}
                gas = CEOSGas(self.class_flash, eos_kwargs=eos_kwargs)
                liquid = CEOSLiquid(self.class_flash, eos_kwargs=eos_kwargs)
                flasher = FlashVL(self.constants, self.properties, liquid=liquid, gas=gas)

                flashed = flasher.flash(T=t[i] + 273.15, P=p[i][j] * pow(10, 5), zs=zs)
                if flashed.VF != 1:  # Есть жидкая фаза
                    for idx, key in enumerate(gas_keys):
                        sum_dict[key] += abs(fractions[key][i][j] - flashed.liquid0.zs[idx]) 
                else:
                    for idx, key in enumerate(gas_keys):
                        sum_dict[key] += abs(fractions[key][i][j] - 0.0)
                count += 1

        avg_error = {key: 0 for key in gas_keys}
        for key in gas_keys:
            avg_error[key] = sum_dict[key] / count

        return sum(avg_error.values())

    def optimize_parameters(self, initial_guess: list[float], fractions: callable, t: list[float],
                            p: list[list[float]], zs: list[float]):
        args = (fractions, t, p, zs)
        bounds = [(-1, 1), (0, 0.004)] * len(zs)
        bounds = tuple(bounds)
        result = minimize(self.func_opt, initial_guess, args=args, bounds=bounds)
        if result.success:
            self.params = result.x
        else:
            raise ValueError(result.message)
        return result

    def calculate_flash(self, t: float, p: float):
        kijs = [[0] * len(self.fractions) for _ in range(len(self.fractions))]
        index = 0
        unique_kijs = set()  # Создаем множество для хранения уникальных значений kijs

        for k in range(len(self.fractions)):
            for l in range(k + 1, len(self.fractions)):
                kijs[k][l] = kijs[l][k] = self.params[index] + (t - self.t_min) * self.params[index + 1]
                unique_kijs.add(kijs[k][l])  # Добавляем значение kijs в множество
                index += 2

        unique_kijs_list = list(unique_kijs)  # Преобразуем множество в список

        eos_kwargs = {'Pcs': self.constants.Pcs, 'Tcs': self.constants.Tcs, 'omegas': self.constants.omegas,
                      'kijs': kijs}
        gas = CEOSGas(self.class_flash, eos_kwargs=eos_kwargs)
        liquid = CEOSLiquid(self.class_flash, eos_kwargs=eos_kwargs)
        flasher = FlashVL(self.constants, self.properties, liquid=liquid, gas=gas)

        PT = flasher.flash(T=t + 273.15, P=p * pow(10, 5), zs=self.fractions)

        if PT.VF == 1:  # Есть только пар
            return [PT.gas.zs, t, p, PT.VF, unique_kijs_list]
        else:  # Не только пар
            return [PT.liquid0.zs, t, p, PT.VF, unique_kijs_list]

'''
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
@dataclass
class Oil:
    M_oil: float
    T_oil: float
    P_oil: float
    Omega_oil: float
    name: str


def create_constants_package(oil: Oil, IDs) -> ChemicalConstantsPackage:
    # Перезаписали критические параметры для воды под масло ISO 32
    constants = ChemicalConstantsPackage.from_IDs(IDs)[0]
    return ChemicalConstantsPackage(
        MWs=constants.MWs[:-1] + [oil.M_oil],
        names=constants.names[:-1] + [oil.name],
        omegas=constants.omegas[:-1] + [oil.Omega_oil],
        Pcs=constants.Pcs[:-1] + [oil.P_oil],
        Tcs=constants.Tcs[:-1] + [oil.T_oil]
    )
