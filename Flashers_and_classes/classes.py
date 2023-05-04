from functions.read_func import read_data_from_file_lines
from thermo import ChemicalConstantsPackage


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
oil_PRMIX = Oil(400.0, 500, 950e3, 0.8)

# Свойства Oil (Молярная масса масла, температура (Кельвины), давление (Па), константа омега)
oil_SRKMIX = Oil(400.0, 500, 950e3, 0.2)

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