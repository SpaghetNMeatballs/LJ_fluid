from dataclasses import dataclass

from density import Density
from functions.chemical_functions import calculate_acentrick_factor
from viscosity import VFT
from common import Point, fill_values
from thermo.chemical import Chemical
from thermo import ViscosityLiquid, VolumeLiquid
from thermo.eos import VDW


@dataclass
class OilTemplate:
    MW: float
    Tc: float
    Pc: float
    omega: float


def generate_oil() -> OilTemplate:
    MW = 400
    Tc = 800
    Pc = 950e3
    omega = calculate_acentrick_factor(Tc, Pc, MW)
    result = OilTemplate(MW=MW, Tc=Tc, Pc=Pc, omega=omega)
    return result


class Material:
    d_method = "LETSOU_STIEL"
    d_method_p = "LUCAS"
    dens_method = "MMSNM0"
    dens_method_p = "COSTALD_COMPRESSED"

    def __init__(
        self,
        inp_chem: Chemical | OilTemplate,
    ):
        self.MW = inp_chem.MW
        self.Tc = inp_chem.Tc
        self.Pc = inp_chem.Pc
        self.Vc = (8.314 * self.Tc) / (8 * self.Pc)
        self.omega = inp_chem.omega
        self.dynamic_visc_calc = ViscosityLiquid(
            MW=self.MW,
            Tc=self.Tc,
            Pc=self.Pc,
            Vc=self.Vc,
            omega=self.omega,
        )
        self.density_calc = VolumeLiquid(
            MW=self.MW,
            Tc=self.Tc,
            Pc=self.Pc,
            Vc=self.Vc,
            omega=self.omega,
        )

    # Возвращает свойства смеси в порядке [динамическая вязкость, кинематическая вязкость, плотность]
    def calc_prop(self, t: float, p: float = None) -> Point:
        t += 273.15
        properties = [None, None, None]
        if p is not None:
            properties[0] = self.dynamic_visc_calc.calculate_P(
                T=t, P=p, method=self.d_method_p
            )
            properties[2] = self.density_calc.calculate_P(
                T=t, P=p, method=self.dens_method_p
            )
        else:
            properties[0] = self.dynamic_visc_calc.calculate(T=t, method=self.d_method)
            properties[2] = self.density_calc.calculate(T=t, method=self.dens_method)
        properties = fill_values(properties)
        result = Point(t=t, p=p, mu=properties[0], v=properties[1], rho=properties[2])
        return result


if __name__ == "__main__":
    r22 = Chemical("r-22")
    oil = generate_oil()
    r22_mat = Material(r22)
    oil_mat = Material(oil)
    print(r22_mat.calc_prop(-10))
