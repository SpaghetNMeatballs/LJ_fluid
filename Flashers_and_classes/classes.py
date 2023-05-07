from functions.read_func import read_data_from_file_lines
from thermo import ChemicalConstantsPackage
from scipy.optimize import minimize
from thermo import CEOSGas, CEOSLiquid, FlashVL


class FlashCalculator:
    def __init__(self, class_flash, properties, constants):
        self.class_flash = class_flash
        self.properties = properties
        self.constants = constants

    def func_opt(self, parametrs, Fraction_R22, Fraction_R115, temperature, pressure, zs):
        count = 0
        sum_R22 = 0
        sum_R115 = 0
        for i in range(len(temperature.values)):
            for j in range(len(pressure.values[i])):
                k13 = parametrs[0] + (temperature.values[i] - temperature.values[0]) * parametrs[1]
                k23 = parametrs[2] + (temperature.values[i] - temperature.values[0]) * parametrs[3]
                k12 = parametrs[4] + (temperature.values[i] - temperature.values[0]) * parametrs[5]
                kijs = [[0, k12, k13],
                        [k12, 0, k23],
                        [k13, k23, 0]]

                eos_kwargs = {'Pcs': self.constants.Pcs, 'Tcs': self.constants.Tcs, 'omegas': self.constants.omegas, 'kijs': kijs}
                gas = CEOSGas(self.class_flash, eos_kwargs=eos_kwargs)
                liquid = CEOSLiquid(self.class_flash, eos_kwargs=eos_kwargs)
                flasher = FlashVL(self.constants, self.properties, liquid=liquid, gas=gas)

                PT = flasher.flash(T=temperature.values[i] + 273.15, P=pressure.values[i][j] * pow(10, 5), zs=zs)
                if PT.VF != 1:
                    ABS_error_R22 = abs(Fraction_R22.values[i][j] - PT.liquid0.zs[0])
                    ABS_error_R115 = abs(Fraction_R115.values[i][j] - PT.liquid0.zs[1])
                else:
                    ABS_error_R22 = abs(Fraction_R22[i][j] - 0.0)
                    ABS_error_R115 = abs(Fraction_R115[i][j] - 0.0)
                count += 1
                sum_R22 += ABS_error_R22
                sum_R115 += ABS_error_R115

        Avg_error_R22 = sum_R22 / count
        Avg_error_R115 = sum_R115 / count
        return Avg_error_R115 + Avg_error_R22

    def optimize_parameters(self, initial_guess_k13_k23, Fraction_R22, Fraction_R115, temperature, pressure, zs):
        result = minimize(self.func_opt, initial_guess_k13_k23, args=(Fraction_R22, Fraction_R115, temperature, pressure, zs), bounds=((-1, 1), (0, 0.004), (-1, 1), (0, 0.004), (-1, 1), (0, 0.004)))
        if result.success:
            print(result.message)
        else:
            raise ValueError(result.message)
        return result

    def calculate_flash(self, T, P, fractions, i, j, result, class_flash, constants, properties):
        k13 = result.x[0] + (T[i] - T[0]) * result.x[1]
        k23 = result.x[2] + (T[i] - T[0]) * result.x[3]
        k12 = result.x[4] + (T[i] - T[0]) * result.x[5]

        kijs = [[0, k12, k13],
                [k12, 0, k23],
                [k13, k23, 0]]

        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
        gas = CEOSGas(class_flash, eos_kwargs=eos_kwargs)
        liquid = CEOSLiquid(class_flash, eos_kwargs=eos_kwargs)
        flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)

        PT = flasher.flash(T=T[i] + 273.15, P=P[i][j] * pow(10, 5), zs=fractions)

        if PT.VF == 1:  # Есть только пар
            return [PT.gas.zs, T, P, PT.VF, k12, k13, k23]
        else:  # Не только пар
            return [PT.liquid0.zs, T, P, PT.VF, k12, k13, k23]


class ChemicalConstants:
    def __init__(self, oil, IDs):
        self.oil = oil
        self.constants, self.properties = ChemicalConstantsPackage.from_IDs(IDs)
        self.constants2 = self.create_constants_package()

    def create_constants_package(self):
        # Перезаписали критические параметры для воды под масло ISO 32
        constants = self.constants
        oil = self.oil
        constants2 = ChemicalConstantsPackage(
            MWs=constants.MWs[:-1] + [oil.M_oil],
            names=constants.names[:-1] + ['ISO 32'],
            omegas=constants.omegas[:-1] + [oil.Omega_oil],
            Pcs=constants.Pcs[:-1] + [oil.P_oil],
            Tcs=constants.Tcs[:-1] + [oil.T_oil]
        )
        return constants2

class Oil:
    def __init__(self, M_oil, T_oil, P_oil, Omega_oil):
        self.M_oil = M_oil
        self.T_oil = T_oil
        self.P_oil = P_oil
        self.Omega_oil = Omega_oil


class PhysicalProperty:
    def __init__(self, values):
        self.values = values

    def __len__(self):
        return len(self.values)

    def __iter__(self):
        return iter(self.values)


class Temperature(PhysicalProperty):
    pass


class Pressure(PhysicalProperty):
    pass


class Viscosity(PhysicalProperty):
    pass


class Density(PhysicalProperty):
    pass


class Fraction(PhysicalProperty):
    pass


# Параметры: массив температур
T = read_data_from_file_lines('../data/temperature.txt')

# Давление R-22: от -10 градусов до 125
pressure_values = read_data_from_file_lines('../data/pressure.txt')

# Вязкость R-22: от -10 градусов до 125
viscosity_values = read_data_from_file_lines('../data/viscosity.txt')

# Плотность R-22:
density_values = read_data_from_file_lines('../data/density.txt')

# Доли R-22 в жидкой фазе
Fraction_R22_values = read_data_from_file_lines('../data/fraction_r22.txt')

# Доли R-115 в жидкой фазе
Fraction_R115_values = read_data_from_file_lines('../data/fraction_r115.txt')

# Свойства Oil (Молярная масса масла, температура (Кельвины), давление (Па), константа омега)
oil_PRMIX = Oil(400.0, 800, 950e3, 0.8)

# Свойства Oil (Молярная масса масла, температура (Кельвины), давление (Па), константа омега)
oil_SRKMIX = Oil(400.0, 800, 950e3, 0.2)

# Создание констант и экземпляров классов
IDs = ['r-22', 'chloropentafluoroethane', 'water']

chemical_constants_PR = ChemicalConstants(oil_PRMIX, IDs).create_constants_package()
properties_PR = ChemicalConstants(oil_PRMIX, IDs)

chemical_constants_SR = ChemicalConstants(oil_SRKMIX, IDs).create_constants_package()
properties_SR = ChemicalConstants(oil_SRKMIX, IDs)

temperature = Temperature(T)
pressure = Pressure(pressure_values)
viscosity = Viscosity(viscosity_values)
density = Density(density_values)
Fraction_R22 = Fraction(Fraction_R22_values)
Fraction_R115 = Fraction(Fraction_R115_values)