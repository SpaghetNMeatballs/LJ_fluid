from thermo import Chemical

from material import Material, generate_oil
from Flashers_and_classes.classes import FlashCalculator
from common import mole_to_mass, aad_new, aad, Point, create_frac_calc_2comp
from data_csv import get_data
from scipy.optimize import minimize
from numpy import array
from viscosity import VFT
from density import Density
from calc_functions import arrhenius, kendall_monroe


class Mixture:
    def __init__(
        self,
        components: list[Material] = None,
        mole_frac_solver: FlashCalculator = None,
        mole_frac: list[float] = None,
        mix_method: callable = None,
    ):
        self.params = []
        self.components = components
        self.mole_frac_solver = mole_frac_solver
        self.mole_frac = mole_frac
        self.mix_method = mix_method
        self.maxiter = 2**20

    def calc_mass_frac(self, point: list[float]) -> list[float]:
        t, p = point
        return mole_to_mass(
            self.mole_frac_solver.calculate_flash(t, p)[0], self.mole_frac
        )

    def calc_all_components_props(self, point: list[float]) -> list[Point]:
        t, p = point
        result = []
        for component in self.components:
            result.append(component.calc_prop(t, p))
        return result

    # Возвращает свойства смеси в порядке
    def calc_mixture_props(self, point: list[float]) -> list[float]:
        t, p = point
        component_props = self.calc_all_components_props([t, p])
        # Группируем свойства в кортежи
        component_props = list(zip(component_props))
        result = [
            self.mix_method(
                self.params[0], [component_props[0][i], component_props[1][i]]
            )
            for i in range(3)
        ]
        for j in range(2, len(component_props)):
            result = [
                self.mix_method(self.params[j - 1], [component_props[j][i], result[i]])
                for i in range(3)
            ]
        return result

    def minimize_params(self, default_params, points: list[list[float]]) -> bool:
        params = [default_params.copy()] * (len(self.components) - 1)
        args = (
            self.calc_mixture_props,
            [[point[0], point[1]] for point in points],
            [point[-1] for point in points],
        )
        result_minimized = minimize(
            fun=aad_new,
            x0=array(params),
            method="Nelder-Mead",
            args=args,
            tol=1e-5,
            options={"maxiter": self.maxiter, "maxfev": self.maxiter},
        )
        return result_minimized.success


if __name__ == "__main__":
    # Loading cavestri data
    print("(!) Loading cavestri data for mixture")
    cavestri_data = get_data(True)
    print("(!) Data loading complete")

    # Mole frac solver
    print("(!) Initialising frac solver")
    mole_frac_solver = create_frac_calc_2comp()
    print("(!) Frac solver initialised")

    # Chemical solvers
    print("(!) Initialising chemical solvers")
    r22_chem = Chemical("r-22")
    r115_chem = Chemical("Chloropentafluoroethane")
    oil_chem = generate_oil()
    print("(!) Chemical solvers initialised")

    # Components
    print("(!) Initialising materials")
    oil = Material(oil_chem)
    r115 = Material(r115_chem)
    r22 = Material(r22_chem)
    print("(!) Materials initialised")

    # Mole fractions for r22 and oil mixture
    g = 0.55
    zs = [g, 1 - g]

    # Mixture
    print("(!) Initialising mixture")
    mixture = Mixture(
        components=[oil, r22], mole_frac_solver=mole_frac_solver, mole_frac=zs
    )
    print("(!) Mixture initialised")
    print("Calculating viscosity in points of cavestri.csv")
    print("T\t\tP\t\tViscExp\t\tViscCalc\tDelta")
    for point in cavestri_data.points:
        t, p, visc = point.t, point.p, point.v
        props = mixture.calc_all_components_props([t, None])
        mole_frac = mixture.mole_frac_solver.calculate_flash(t, point.p)[0]
        for i in range(len(props)):
            props[i].x = mole_frac[i]
        visc_calc = kendall_monroe(props)
        print(
            f"{t:.2f}\t{p:.2f}\t{visc:.2f}\t\t{visc_calc:.2f}\t\t{aad([visc], [visc_calc]):.2f}%"
        )
