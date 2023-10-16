from material import Material
from Flashers_and_classes.classes import FlashCalculator
from common import mole_to_mass


class Mixture:
    def __init__(self, components: list[Material] = None, mole_frac_solver: FlashCalculator = None,
                 mole_frac: list[float] = None,
                 mix_method: callable[[list[float], list[float], list[float]], float] = None):
        self.params = []
        self.components = components
        self.mole_frac_solver = mole_frac_solver
        self.mole_frac = mole_frac
        self.mix_method = mix_method

    def calc_mass_frac(self, t: float, p: float) -> list[float]:
        return mole_to_mass(self.mole_frac_solver.calculate_flash(t, p), self.mole_frac)

    def calc_all_components_props(self, t: float, p: float = None) -> list[list[float]]:
        result = []
        for component in self.components:
            result.append(component.calc_prop(t, p))
        return result

    # Возвращает свойства смеси в порядке
    def calc_mixture_props(self, t: float, p: float = None) -> list[float]:
