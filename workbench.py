from typing import List
from scipy.optimize import minimize
import pandas as pd
from calc_functions import arrhenius
from data_csv import DataManager, get_data
from viscosity import VFT
import numpy as np


def load_mole_fractions_from_csv(path='flasher_data/flash_results_dataset.csv'):
    df = pd.read_csv(path)
    return df[['T_C', 'P_bar', 'x_R22', 'x_R115', 'x_ISO32']]


def merge_data(visc_data: DataManager, mole_frac_df: pd.DataFrame) -> DataManager:
    df_visc = pd.DataFrame({
        'T': visc_data.T,
        'P': visc_data.P,
        'visc': visc_data.visc,
        'rho': visc_data.rho
    })

    df_visc = df_visc.sort_values(['T', 'P']).reset_index(drop=True)
    mole_frac_df_sorted = mole_frac_df.sort_values(['T_C', 'P_bar']).reset_index(drop=True)

    merged = pd.merge_asof(df_visc, mole_frac_df_sorted, left_on='T', right_on='T_C', direction='nearest', tolerance=1.0)

    combined_data = DataManager()
    for _, row in merged.iterrows():
        if pd.isna(row['x_R22']) or pd.isna(row['x_R115']):
            continue

        liq_frac1 = 0.0
        liq_frac2 = 0.0

        combined_data.add_point([
            row['T'],
            row['P'],
            row['x_R22'],
            row['x_R115'],
            liq_frac1,
            liq_frac2,
            row['visc'],
            row['rho']
        ])
    return combined_data


def loss_arrhenius(params, T_list, mix_visc_list, fractions, vft_r22, vft_r115):
    errors = []
    for T, exp_visc, (frac_r22, frac_r115) in zip(T_list, mix_visc_list, fractions):
        mu1 = vft_r22.calc_visc_def(T)
        mu2 = vft_r115.calc_visc_def(T)
        pred = arrhenius([mu1, mu2], params)
        errors.append((pred - exp_visc) ** 2)
    return sum(errors) / len(errors)


def load_flash_csv_data(path='flasher_data/flash_results_dataset.csv') -> DataManager:
    df = pd.read_csv(path)
    data = DataManager()
    for idx, row in df.iterrows():
        T = float(row['T_C'])
        P = float(row['P_bar'])
        x_R22 = float(row['x_R22'])
        x_R115 = float(row['x_R115'])
        liq_frac1 = 0.0
        liq_frac2 = 0.0
        visc = 0.0
        rho = 0.0
        data.add_point([T, P, x_R22, x_R115, liq_frac1, liq_frac2, visc, rho])
    return data


def loss_arrhenius(params: List[float], T_list: List[float], mix_visc_list: List[float], fractions: List[List[float]],
                   vft_r22: VFT, vft_r115: VFT) -> float:
    errors = []
    for T, exp_visc, (frac_r22, frac_r115) in zip(T_list, mix_visc_list, fractions):
        mu1 = vft_r22.calc_visc_def(T)
        mu2 = vft_r115.calc_visc_def(T)
        pred = arrhenius([mu1, mu2], params)
        if exp_visc > 0:
            errors.append((np.log(pred) - np.log(exp_visc)) ** 2)
    return sum(errors) / len(errors) if errors else 1e6


def build_mu_data(data: DataManager, component_index: int):
    temps = []
    viscs = []
    for T, fracs, mu in zip(data.T, data.gas_fractions, data.visc):
        x = fracs[component_index]
        if x == 0:
            continue
        temps.append(T)
        viscs.append(mu / x)
    return temps, viscs


def prepare_vft_and_fit(data: DataManager):
    temps_r22, viscs_r22 = build_mu_data(data, component_index=0)
    temps_r115, viscs_r115 = build_mu_data(data, component_index=1)

    x0 = [1e-5, 1000, 200]

    vft_r22 = VFT(temps_r22, viscs_r22, x0.copy())
    vft_r115 = VFT(temps_r115, viscs_r115, x0.copy())

    initial_guess = [0.5, 0.5]
    bounds = [(0, 2), (0, 2)]
    constraints = {'type': 'eq', 'fun': lambda x: x[0] + x[1] - 1}

    res = minimize(loss_arrhenius, initial_guess, args=(data.T, data.visc, data.gas_fractions, vft_r22, vft_r115),
                   bounds=bounds, constraints=constraints,
                   options={'maxfev': 10000, 'disp': True})

    if not res.success:
        print("Optimization failed:", res.message)
        x1, x2 = 0.5, 0.5
    else:
        x1, x2 = res.x

    return vft_r22, vft_r115, x1, x2


def print_results(data: DataManager, vft_r22: VFT, vft_r115: VFT, x1: float, x2: float):
    print(f"\nОптимальные параметры Аррениуса: x1 = {x1:.4f}, x2 = {x2:.4f}\n")
    print(
        f"{'T':>6} {'P':>6} {'mu_R22':>10} {'mu_R115':>10} {'rho_R22':>10} {'rho_R115':>10} "
        f"{'mu_mix(pred)':>14} {'mu_mix(exp)':>12} {'rho_mix':>10} {'Abs_err':>10} {'Rel_err,%':>10}"
    )
    print("-" * 110)

    for T, P, fracs, rho, mu_exp in zip(data.T, data.P, data.gas_fractions, data.rho, data.visc):
        x1_gas, x2_gas = fracs
        if x1_gas == 0 or x2_gas == 0:
            continue

        mu_r22 = vft_r22.calc_visc_def(T)
        mu_r115 = vft_r115.calc_visc_def(T)
        mu_mix = arrhenius([mu_r22, mu_r115], [x1, x2])

        rho_r22 = rho / x1_gas if x1_gas else 0
        rho_r115 = rho / x2_gas if x2_gas else 0

        abs_err = abs(mu_mix - mu_exp)
        rel_err = (abs_err / mu_exp * 100) if mu_exp != 0 else 0.0

        print(
            f"{T:6.1f} {P:6.2f} {mu_r22:10.4f} {mu_r115:10.4f} {rho_r22:10.4f} {rho_r115:10.4f} "
            f"{mu_mix:14.4f} {mu_exp:12.4f} {rho:10.4f} {abs_err:10.4f} {rel_err:10.2f}"
        )

def data_manager_pipeline():
    data = get_data()
    vft_r22, vft_r115, x1, x2 = prepare_vft_and_fit(data)
    print_results(data, vft_r22, vft_r115, x1, x2)


def flasher_csv_pipeline():
    data = load_flash_csv_data()
    vft_r22, vft_r115, x1, x2 = prepare_vft_and_fit(data)
    print_results(data, vft_r22, vft_r115, x1, x2)


def combined_pipeline():
    visc_data = get_data()
    mole_frac_df = load_mole_fractions_from_csv()

    combined_data = merge_data(visc_data, mole_frac_df)

    vft_r22, vft_r115, x1, x2 = prepare_vft_and_fit(combined_data)

    print_results(combined_data, vft_r22, vft_r115, x1, x2)


if __name__ == "__main__":
    combined_pipeline()
