from material import Material
from pytest import fixture
from density import Density
from viscosity import VFT


@fixture
def no_solvers():
    return Material(assign_flag=False)


@fixture
def v_solver():
    return Material(assign_flag=True, d_visc_calc=VFT())


@fixture
def k_solver():
    return Material(assign_flag=True, k_visc_calc=VFT())


@fixture
def dens_solver():
    return Material(assign_flag=True, dens_calc=Density(0, 0, 0))


@fixture
def dens_v_solver():
    return Material(assign_flag=True, dens_calc=Density(0, 0, 0), d_visc_calc=VFT())


@fixture
def dens_k_solver():
    return Material(assign_flag=True, dens_calc=Density(0, 0, 0), k_visc_calc=VFT())


@fixture
def full_solver():
    return Material(assign_flag=True, k_visc_calc=VFT(), d_visc_calc=VFT(), dens_calc=Density(0, 0, 0))


def test_calc_amount(no_solvers, dens_solver, k_solver, v_solver, dens_k_solver, dens_v_solver, full_solver):
    assert no_solvers.get_calc_amount() == 0
    assert dens_solver.get_calc_amount() == 1
    assert k_solver.get_calc_amount() == 1
    assert v_solver.get_calc_amount() == 1
    assert dens_k_solver.get_calc_amount() == 2
    assert dens_v_solver.get_calc_amount() == 2
    assert full_solver.get_calc_amount() == 3
