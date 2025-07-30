import matplotlib.pyplot as plt
from thermo import PRMIX
from thermo.chemical import Chemical
from Flashers_and_classes.classes import FlashCalculator, create_constants_package, Oil
from functions.chemical_functions import mass_fractions_to_mole_fractions, calculate_acentrick_factor
from functions.read_func import read_data_from_file_lines
import openpyxl
import csv


def get_chemical_molar_masses(components):
    return {name: round(Chemical(name).MW, 2) for name in components}


def prepare_experimental_data():
    temperature = read_data_from_file_lines('temperature.txt')
    pressure = read_data_from_file_lines('pressure.txt')
    r22_fractions = read_data_from_file_lines('fraction_r22.txt')
    r115_fractions = read_data_from_file_lines('fraction_r115.txt')

    return temperature, pressure, {
        'R-22': r22_fractions,
        'R-115': r115_fractions
    }


def compute_mole_fractions(fractions, molar_masses, molar_mass_oil):
    return mass_fractions_to_mole_fractions(fractions, molar_masses, molar_mass_oil)


def build_flash_calculator(zs, t_min, molar_mass_oil, omega, ids):
    oil = Oil(molar_mass_oil, 800, 950e3, omega, 'ISO32')
    constants = create_constants_package(oil, ids)
    initial_k = [0.21, 0.001, 0.35, 0.001, 0.1, 0.001]
    return FlashCalculator(PRMIX, oil, constants, initial_k, zs, t_min), initial_k


def write_headers(sheet):
    sheet.append([
        'Т(°C)', 'P(кПа)', 'x_R-22(модел.)', 'х_R-115(модел.)', 'х_ISO 32(модел.)',
        'x_R-22(эксп.)', 'x_R-115(эксп.)', 'x_ISO 32(эксп.)',
        'Aбс.ошибка R-22', 'Абс.ошибка R-115'
    ])


def process_flash_results(flash_calculator, T, P, mole_fractions_exp, zs_mole, sheet):
    total_points, total_error_r22, total_error_r115 = 0, 0, 0

    for i, temp in enumerate(T):
        r22_model, r115_model, r22_exp, r115_exp, p_plot = [], [], [], [], []

        error_r22_T, error_r115_T = 0, 0

        for j, pressure in enumerate(P[i]):
            flash_data = flash_calculator.calculate_flash(temp, pressure)

            if flash_data[3] == 1:
                flash_data[0][0] = 0.0
                flash_data[0][1] = 0.0

            abs_err_r22 = abs(mole_fractions_exp['R-22'][i][j] - flash_data[0][0])
            abs_err_r115 = abs(mole_fractions_exp['R-115'][i][j] - flash_data[0][1])

            sheet.append([
                temp, pressure, flash_data[0][0], flash_data[0][1], flash_data[0][2],
                mole_fractions_exp['R-22'][i][j], mole_fractions_exp['R-115'][i][j],
                mole_fractions_exp['Oil'][i][j], abs_err_r22, abs_err_r115
            ])

            r22_model.append(flash_data[0][0])
            r115_model.append(flash_data[0][1])
            r22_exp.append(mole_fractions_exp['R-22'][i][j])
            r115_exp.append(mole_fractions_exp['R-115'][i][j])
            p_plot.append(pressure)

            error_r22_T += abs_err_r22
            error_r115_T += abs_err_r115
            total_error_r22 += abs_err_r22
            total_error_r115 += abs_err_r115
            total_points += 1

        plot_results(temp, p_plot, r22_model, r115_model, r22_exp, r115_exp)

    avg_r22 = total_error_r22 / total_points
    avg_r115 = total_error_r115 / total_points

    print(f'\n****** Средняя погрешность ******')
    print(f'R-22: {avg_r22:.5f} ({avg_r22 * 100:.2f}%)')
    print(f'R-115: {avg_r115:.5f} ({avg_r115 * 100:.2f}%)')

    return sheet


def plot_results(temp, pressures, r22_model, r115_model, r22_exp, r115_exp):
    plt.title(f'Зависимость x от P при T={temp}°C', fontsize=16)
    plt.xlabel('P, бар', fontsize=14)
    plt.ylabel('x, мольн. доля', fontsize=14)

    plt.scatter(pressures, r22_model)
    plt.scatter(pressures, r115_model)
    plt.scatter(pressures, r22_exp)
    plt.scatter(pressures, r115_exp)

    plt.plot(pressures, r22_model)
    plt.plot(pressures, r115_model)

    plt.legend(['R-22 модел.', 'R-115 модел.', 'R-22 эксп.', 'R-115 эксп.'])
    plt.grid(True)
    plt.show()


def export_to_csv(sheet, filepath):
    with open(filepath, mode='w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
        writer.writerow([
            'T_C', 'P_bar', 'x_R22', 'x_R115', 'x_ISO32',
            'x_R22_exp', 'x_R115_exp', 'x_ISO32_exp', 'ABS_err_R22', 'ABS_err_R115'
        ])
        for row in sheet.iter_rows(min_row=2, values_only=True):
            writer.writerow(row)


def main():
    # Исходные данные
    gas_mass_ratio = 0.7
    zs_mass = [0.488 * gas_mass_ratio, 0.512 * gas_mass_ratio, 1 - gas_mass_ratio]
    ids = ['r-22', 'chloropentafluoroethane', 'water']
    M_oil, T_oil, P_oil = 400.0, 800, 950e3

    # Чтение и преобразование данных
    T, P, exp_mass_fractions = prepare_experimental_data()
    t_min = T[0]

    molar_masses = get_chemical_molar_masses(ids[:2])
    molar_masses_dict = {'R-22': molar_masses[ids[0]], 'R-115': molar_masses[ids[1]]}
    mole_fractions_exp = compute_mole_fractions(exp_mass_fractions, molar_masses_dict, M_oil)

    zs_mass_dict = {
        'R-22': [[zs_mass[0]]],
        'R-115': [[zs_mass[1]]]
    }
    zs_mole = compute_mole_fractions(zs_mass_dict, molar_masses_dict, M_oil)
    zs_mole_values = [zs_mole['R-22'][0][0], zs_mole['R-115'][0][0], zs_mole['Oil'][0][0]]

    omega = calculate_acentrick_factor(T_oil, P_oil, M_oil)
    flash_calculator, initial_guess = build_flash_calculator(zs_mole_values, t_min, M_oil, omega, ids)
    flash_calculator.optimize_parameters(initial_guess, mole_fractions_exp, T, P, zs_mole_values)

    workbook = openpyxl.Workbook()
    sheet = workbook.active
    write_headers(sheet)

    process_flash_results(flash_calculator, T, P, mole_fractions_exp, zs_mole_values, sheet)
    export_to_csv(sheet, "..\\flasher_data\\flash_results_dataset.csv")


if __name__ == "__main__":
    main()
