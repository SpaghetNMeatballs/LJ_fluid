from material import Material
from Flashers_and_classes.classes import FlashCalculator
from common import mole_to_mass, aad_new, aad, Point, create_frac_calc_2comp
from data_csv import get_data
from scipy.optimize import minimize
from numpy import array
from viscosity import VFT
from density import Density
from calc_functions import arrhenius


class Mixture:
    def __init__(self, components: list[Material] = None, mole_frac_solver: FlashCalculator = None,
                 mole_frac: list[float] = None, mix_method: callable = None):
        self.params = []
        self.components = components
        self.mole_frac_solver = mole_frac_solver
        self.mole_frac = mole_frac
        self.mix_method = mix_method
        self.maxiter = 2 ** 20

    def calc_mass_frac(self, point: list[float]) -> list[float]:
        t, p = point
        return mole_to_mass(self.mole_frac_solver.calculate_flash(t, p)[0], self.mole_frac)

    def calc_all_components_props(self, point: list[float]) -> list[list[float]]:
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
        result = [self.mix_method(self.params[0], [component_props[0][i], component_props[1][i]]) for i in range(3)]
        for j in range(2, len(component_props)):
            result = [self.mix_method(self.params[j - 1], [component_props[j][i], result[i]]) for i in range(3)]
        return result

    def minimize_params(self, default_params, points: list[list[float]]) -> bool:
        params = [default_params.copy()] * (len(self.components) - 1)
        args = (self.calc_mixture_props, [[point[0], point[1]] for point in points], [point[-1] for point in points])
        result_minimized = minimize(fun=aad_new, x0=array(params), method='Nelder-Mead', args=args, tol=1e-5,
                                    options={'maxiter': self.maxiter, 'maxfev': self.maxiter})
        return result_minimized.success


if __name__ == '__main__':

    # Loading cavestri data
    print('(!) Loading cavestri data for mixture')
    cavestri_data = get_data(True)
    print('(!) Data loading complete')

    # Mole frac solver
    print('(!) Initialising frac solver')
    mole_frac_solver = create_frac_calc_2comp()
    print('(!) Frac solver initialised')

    # Viscosity solvers for oil, r115 and r22
    print('(!) Initialising visc solvers')
    vft_oil = VFT()
    vft_oil.params = [0.005337316178523319, 1755.4979644152024, 109.82206931834386, 0, 0, 0, 0]
    vft_r115 = VFT()
    vft_r115.params = [0.0008577907640623858, -69.08004012663991, 4.2625169264404175, 0, 0, 0, 0]
    vft_r22 = VFT()
    vft_r22.params = [2.2755471578933423e-05, 626.8283963685842, 2.50475798005479, 0, 0, 0, 0]
    print('(!) Visc solvers initialised')

    # P0 for density solvers
    point0_r22 = Point(t=-60, p=0.07, rho=1463.67431328749)
    point0_r115 = Point(t=-49.5, p=0.07, rho=1581.88691446926)
    point0_oil = Point(t=-25, rho=0.935)

    # Density solvers
    print('(!) Initialising density solvers')
    density_r22 = Density(point0_r22)
    density_r115 = Density(point0_r115)
    density_oil = Density(point0_oil)
    print('(!) Density solvers initialised')

    # Components
    print('(!) Initialising materials')
    oil = Material(d_visc_calc=vft_oil, dens_calc=density_oil)
    r115 = Material(d_visc_calc=vft_r115, dens_calc=density_r115)
    r22 = Material(d_visc_calc=vft_r22, dens_calc=density_r22)
    print('(!) Materials initialised')

    # Mole fractions for r22 and oil mixture
    g = 0.55
    zs = [g, 1 - g]

    # Mixture
    print('(!) Initialising mixture')
    mixture = Mixture(components=[oil, r22], mole_frac_solver=mole_frac_solver, mole_frac=zs)
    print('(!) Mixture initialised')
    print('Calculating viscosity in points of cavestri.csv')
    print('T\t\tP\t\tViscExp\t\tViscCalc\tDelta')
    for point in cavestri_data.points:
        t, p, visc = point.t, point.p, point.v
        props = mixture.calc_all_components_props([t, None])
        mole_frac = mixture.mole_frac_solver.calculate_flash(t, point.p)[0]
        visc_calc = arrhenius([props[0][0], props[1][0], mole_frac[1], mole_frac[0]])
        print(f'{t:.2f}\t{p:.2f}\t{visc:.2f}\t\t{visc_calc:.2f}\t\t{aad([visc], [visc_calc]):.2f}%')
