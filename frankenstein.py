import matplotlib.pyplot as plt
from numpy import arange, array, zeros
from scipy.optimize import minimize
import matplotlib
from thermo import Chemical

from calc_functions import lederer_roegiers, grunberg_nissan, arrhenius, kendall_monroe
from common import Point, create_frac_calc_2comp, mole_to_mass, aad
from data_csv import get_data
from density import Density
from mixture import Mixture
from material import Material, generate_oil
from viscosity import VFT


def mix_visc_check():
    t_list = [-10, 0, 20, 40, 70, 80, 125]
    display_t = 60
    step = 0.1
    """visc_r115 = data_r115.visc_d
    visc_r22 = data_r22.visc_d

    visc_oil = {pref: {-15.171: 413.005, -10.1111: 277.596, -5.0512: 191.427, 0.00868263: 137.18, 4.87395: 99.5737,
                       9.93383: 75.1096, 14.9937: 57.3869, 20.0536: 44.4115, 24.9189: 35.2622, 29.9787: 28.359,
                       35.0386: 23.1013, 39.9039: 19.0611, 44.9638: 15.9304, 50.0237: 13.4856, 54.8889: 11.5632,
                       59.9488: 9.91492, 65.0087: 8.7223, 70.0686: 7.67312, 74.9338: 6.75015, 79.9937: 6.01479,
                       85.0536: 5.35955, 90.1135: 4.83728, 95.1734: 4.31031, 100.039: 3.94046}}"""
    # Viscosity solvers for oil, r115 and r22
    vft_oil = VFT()
    vft_oil.params = [
        0.005337316178523319,
        1755.4979644152024,
        109.82206931834386,
        0,
        0,
        0,
        0,
    ]
    vft_r115 = VFT()
    vft_r115.params = [
        0.0008577907640623858,
        -69.08004012663991,
        4.2625169264404175,
        0,
        0,
        0,
        0,
    ]
    vft_r22 = VFT()
    vft_r22.params = [
        2.2755471578933423e-05,
        626.8283963685842,
        2.50475798005479,
        0,
        0,
        0,
        0,
    ]

    point0_r22 = Point(t=-60, p=0.07, rho=1463.67431328749)
    point0_r115 = Point(t=-49.5, p=0.07, rho=1581.88691446926)
    point0_oil = Point(t=-25, rho=0.935)

    density_r22 = Density(point0_r22)
    density_r115 = Density(point0_r115)
    density_oil = Density(point0_oil)

    # Присваиваем параматреы, расчитанные ранее
    density_r22.e, density_r22.beta = 132772382183158.875000, 0.001752
    density_r115.e, density_r115.beta = 21106748519192.570312, 0.001916
    density_oil.beta = 0.0010332054833384788

    print("Initializing mole frac solver")
    mole_frac_solver = create_frac_calc_2comp()
    print("Mole frac solver setup complete")

    print("Loading cavestri data for mixture")
    cavestri_data = get_data(True)
    print("Data loading complete")
    points_mu_t = cavestri_data.get_by_t(display_t)
    p_scatter = [i.p for i in points_mu_t]
    mu_scatter = [i.mu for i in points_mu_t]

    r22_chem = Chemical("r-22")
    r115_chem = Chemical("Chloropentafluoroethane")
    oil_chem = generate_oil()
    oil = Material(oil_chem)
    r115 = Material(r115_chem)
    r22 = Material(r22_chem)
    g = 0.55
    zs = [g, 1 - g]
    mixture = Mixture(
        components=[oil, r22], mole_frac_solver=mole_frac_solver, mole_frac=zs
    )
    thermo_points = []

    print("Initializing calculation of viscosity points")

    p_calc = arange(min(p_scatter), max(p_scatter) + step, step)
    visc_calc = []
    points_minimize = []
    exp_minimize = []
    for i in cavestri_data.points:
        mole_frac = mole_frac_solver.calculate_flash(i.t, i.p)[0]

        mu_r22 = vft_r22.calc_visc_def(i.t)

        v_oil_p = vft_oil.calc_visc_def(i.t)
        rho_oil_p = density_oil.calc_rho_t(i.t)
        mu_oil = v_oil_p * rho_oil_p

        points_minimize.append([mu_r22, mu_oil, mole_frac[0], mole_frac[1]])
        exp_minimize.append(i.mu)

    args_lr = (lederer_roegiers, points_minimize, exp_minimize)
    args_gn = (grunberg_nissan, points_minimize, exp_minimize)

    for p in p_calc:
        cur_point = Point(t=display_t, p=p * 10)
        mole_frac = mole_frac_solver.calculate_flash(display_t, p)[0]
        mass_frac = mole_to_mass(mole_frac, [86.47, 400.0])

        mu_r22_p = vft_r22.calc_visc_def(display_t)
        rho_r22_p = density_r22.calc_rho_full(cur_point)
        v_r22_p = mu_r22_p / rho_r22_p * 10**6
        point_r22 = Point(
            t=cur_point.t,
            p=cur_point.p,
            mu=mu_r22_p,
            rho=rho_r22_p,
            v=v_r22_p,
            x=mole_frac[0],
        )

        """mu_r115_p = vft_r115.calc_visc_def(display_t)
        rho_r115_p = density_r115.calc_rho_full(cur_point)
        v_r115_p = mu_r115_p / rho_r115_p"""

        v_oil_p = vft_oil.calc_visc_def(display_t)
        rho_oil_p = density_oil.calc_rho_t(display_t)
        mu_oil_p = v_oil_p * rho_oil_p
        point_oil = Point(
            t=cur_point.t,
            p=cur_point.p,
            mu=mu_oil_p,
            rho=rho_oil_p,
            v=v_oil_p,
            x=mole_frac[1],
        )

        # mu_refr_list = [mu_r22_p, mu_r115_p]
        # mu_refr = arrhenius(lr_params[0], mu_refr_list, normalise(mole_frac[:2]))
        visc_calc.append(kendall_monroe([point_r22, point_oil]))
        thermo_props = mixture.calc_all_components_props([cur_point.t, None])
        for i in range(len(thermo_props)):
            thermo_props[i].x = mole_frac[i]
        thermo_points.append(kendall_monroe(thermo_props))
        # visc_calc.append(lederer_roegiers(alpha, [mu_r22_p, mu_oil_p, mole_frac[0], mole_frac[1]]))
    print("Calculation complete, initiating plotting")

    for t_cur in t_list:
        t_cur_points = cavestri_data.get_by_t(t_cur)
        if len(t_cur_points) == 0:
            continue
        points = [i.p for i in t_cur_points]
        x_exp = [i.mu for i in t_cur_points]
        x_calc = []
        for i in t_cur_points:
            mole_frac = mole_frac_solver.calculate_flash(t_cur, i.p)[0]
            mass_frac = mole_to_mass(mole_frac, [86.47, 400.0])

            mu_r22_p = vft_r22.calc_visc_def(t_cur)
            rho_r22_p = density_r22.calc_rho_full(i)
            v_r22_p = mu_r22_p / rho_r22_p
            point_r22 = Point(
                t=i.t,
                p=i.p,
                mu=mu_r22_p,
                rho=rho_r22_p,
                v=v_r22_p,
                x=mole_frac[0],
            )

            # mu_r115_p = vft_r115.calc_visc_def(display_t)
            # rho_r115_p = density_r115.calc_rho_full(cur_point)
            # v_r115_p = mu_r115_p / rho_r115_p

            v_oil_p = vft_oil.calc_visc_def(t_cur)
            rho_oil_p = density_oil.calc_rho_t(t_cur)
            mu_oil_p = v_oil_p * rho_oil_p
            point_oil = Point(
                t=i.t, p=i.p, mu=mu_oil_p, rho=rho_oil_p, v=v_oil_p, x=mole_frac[1]
            )

            # mu_refr_list = [mu_r22_p, mu_r115_p]
            # mu_refr = arrhenius(lr_params[0], mu_refr_list, normalise(mole_frac[:2]))
            # visc_calc.append(arrhenius(lr_params[1], [mu_refr, mu_oil_p], [sum(mole_frac[:2]), mole_frac[2]]))
            x_calc.append(arrhenius([point_r22, point_oil]))
        print(f"AAD для t={t_cur}:\t{aad(x_exp, x_calc)}")

    font = {"family": "Times New Roman", "size": 20}
    matplotlib.rc("font", **font)
    fig, ax = plt.subplots()
    ax.scatter(p_scatter, mu_scatter, label="Экспериментальные данные")
    ax.plot(p_calc, visc_calc, "r", label="Расчётные значения")
    ax.plot(p_calc, thermo_points, "g", label="Значения thermo")
    # ax.plot(arange(temps[0], temps[-1] + 1), viscs_calc, "r",
    #        label="Расчётные значения")
    ax.set_ylabel("$\it{Вязкость, сП}$")
    ax.set_xlabel("$\it{P, бар}$")
    plt.grid()
    plt.legend()
    plt.show()


if __name__ == "__main__":
    mix_visc_check()
