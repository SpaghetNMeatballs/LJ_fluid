import csv
from math import log, exp
from typing import List
from scipy.optimize import minimize
import numpy as np
import pandas as pd
from calc_functions import arrhenius
from data_csv import DataManager, get_data
from viscosity import VFT

def loss_arrhenius(params, T_list, mix_visc_list, fractions, vft_r22, vft_r115):
    errors = []
    for T, exp_visc, (frac_r22, frac_r115) in zip(T_list, mix_visc_list, fractions):
        mu1 = vft_r22.calc_visc_def(T)
        mu2 = vft_r115.calc_visc_def(T)
        pred = arrhenius([mu1, mu2], params)
        errors.append((pred - exp_visc) ** 2)
    return sum(errors) / len(errors)

def load_flash_csv_data(path='Flashers_and_classes\\flash_results_dataset.csv'):
    df = pd.read_csv(path)

    class DataContainer:
        def __init__(self, df):
            self.T = df['T_C'].tolist()
            self.P = df['P_bar'].tolist()
            self.visc = df['mu_cP'].tolist()
            self.rho = df['rho_g_per_ml'].tolist()
            self.gas_fractions = list(zip(df['x_R22'], df['x_R115']))

    return DataContainer(df)

def fit_arrhenius_params(data: DataManager, vft_r22: VFT, vft_r115: VFT):
    T_list = data.T
    mix_visc_list = data.visc
    fractions = data.gas_fractions

    initial_guess = [0.5, 0.5]
    bounds = [(0, 1), (0, 1)]
    constraints = {'type': 'eq', 'fun': lambda x: x[0] + x[1] - 1}

    res = minimize(loss_arrhenius, initial_guess, args=(T_list, mix_visc_list, fractions, vft_r22, vft_r115),
                   bounds=bounds, constraints=constraints)

    return res.x

def build_mu_data(data: DataManager, component_index: int) -> dict:
    mu_data = {}
    for T, P, fracs, mu in zip(data.T, data.P, data.gas_fractions, data.visc):
        x = fracs[component_index]
        if x == 0:
            continue
        mu_individual = mu / x
        if P not in mu_data:
            mu_data[P] = {}
        mu_data[P][T] = mu_individual
    return mu_data

def prepare_vft_and_fit(data):
    mu_data_r22 = build_mu_data(data, component_index=0)
    mu_data_r115 = build_mu_data(data, component_index=1)

    pref_pressure = min(mu_data_r22.keys(), key=lambda p: abs(p - 0.62))
    vft_r22 = VFT(mu_data_r22, pref=pref_pressure)

    pref_pressure_115 = min(mu_data_r115.keys(), key=lambda p: abs(p - 0.62))
    vft_r115 = VFT(mu_data_r115, pref=pref_pressure_115)

    assert vft_r22.setup_def(), "Ошибка подбора VFT для R22"
    assert vft_r115.setup_def(), "Ошибка подбора VFT для R115"

    x1, x2 = fit_arrhenius_params(data, vft_r22, vft_r115)
    return vft_r22, vft_r115, x1, x2


def print_results(data, vft_r22, vft_r115, x1, x2):
    print(f"\nОптимальные параметры Аррениуса: x1 = {x1:.4f}, x2 = {x2:.4f}\n")
    print(f"{'T':>6} {'P':>6} {'mu_R22':>10} {'mu_R115':>10} {'rho_R22':>10} {'rho_R115':>10} {'mu_mix(pred)':>14} {'mu_mix(exp)':>12} {'rho_mix':>10}")
    print("-" * 90)

    for T, P, fracs, rho, mu_exp in zip(data.T, data.P, data.gas_fractions, data.rho, data.visc):
        x1_gas, x2_gas = fracs
        if x1_gas == 0 or x2_gas == 0:
            continue

        mu_r22 = vft_r22.calc_visc_def(T)
        mu_r115 = vft_r115.calc_visc_def(T)
        mu_mix = arrhenius([mu_r22, mu_r115], [x1, x2])

        rho_r22 = rho * (1 / x1_gas)
        rho_r115 = rho * (1 / x2_gas)

        print(f"{T:6.1f} {P:6.2f} {mu_r22:10.4f} {mu_r115:10.4f} {rho_r22:10.4f} {rho_r115:10.4f} {mu_mix:14.4f} {mu_exp:12.4f} {rho:10.4f}")


def data_manager_pipeline():
    data = get_data()
    vft_r22, vft_r115, x1, x2 = prepare_vft_and_fit(data)
    print_results(data, vft_r22, vft_r115, x1, x2)


def flasher_csv_pipeline():
    data = load_flash_csv_data()
    vft_r22, vft_r115, x1, x2 = prepare_vft_and_fit(data)
    print_results(data, vft_r22, vft_r115, x1, x2)



if __name__ == "__main__":
    flasher_csv_pipeline()