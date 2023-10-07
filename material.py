from density import Density
from viscosity import VFT
from common import Point


class Material:
    def get_calc_amount(self) -> int:
        result = 0
        if self.d_viscosity is not None:
            result += 1
        if self.k_viscosity is not None:
            result += 1
        if self.density is not None:
            result += 1
        return result

    def __assign_calculators(self, d_visc_calc: VFT = None, k_visc_calc: VFT = None, dens_calc: Density = None):
        if d_visc_calc:
            self.d_viscosity = d_visc_calc
        if k_visc_calc:
            self.k_viscosity = k_visc_calc
        if dens_calc:
            self.density = dens_calc

    def __init__(self, assign_flag: bool = False, d_visc_calc: VFT = None, k_visc_calc: VFT = None,
                 dens_calc: Density = None):
        self.d_viscosity = None
        self.k_viscosity = None
        self.density = None
        if assign_flag:
            self.__assign_calculators(d_visc_calc, k_visc_calc, dens_calc)

    # Fills missing values of [d_visc, k_visc, dens] list
    def fill_values(self, val_list: list[float | None]) -> list[float]:
        if None not in val_list:
            return val_list
        none_ind = val_list.index(None)
        match none_ind:
            case 0:
                val_list[0] = val_list[1] * val_list[2]
            case 1:
                val_list[1] = val_list[0] / val_list[2]
            case 2:
                val_list[2] = val_list[]


    # Returns list of [d_visc, k_visc, density]
    def calc_prop(self, t: float, p: float = None) -> list[float]:
        if self.get_calc_amount() < 2:
            raise Exception('Not enought solvers initialised')
        result = [None, None, None]
        if self.d_viscosity is not None:
            if p is not None:
                result[0] = self.d_viscosity.calc_visc(t, p)
            else:
                result[0] = self.d_viscosity.calc_visc_def(t)
        if self.k_viscosity is not None:
            if p is not None:
                result[1] = self.k_viscosity.calc_visc(t, p)
            else:
                result[1] = self.k_viscosity.calc_visc_def(t)
        if self.density is not None:
            result[2] = self.density.calc_rho_full(Point(t=t, p=p))

        transfer_funcs = [lambda x: x[1] * x[2], lambda x: x[0] / x[2], lambda x: x[1] / x[0]]
        if None in result:
            none_ind = result.index(None)
            result[none_ind] = transfer_funcs[none_ind](result)
        return result
