# Задача. Предсказывание состава жидкой фазы смеси.
# Входные параметры: Температура C, давление Бар (МПа/10) для каждой точки, вязкость cP, плотность г/мл и состав газа
# R-22/R-115 (48.8/51.2) ------> (0.488/0.512)
# Построить зависимость концентрации от давления

import matplotlib.pyplot as plt
from thermo import PRMIX
from Flashers_and_classes.classes import FlashCalculator, create_constants_package, Oil
from functions.chemical_functions import mass_fractions_to_mole_fractions, calculate_acentrick_factor
from functions.read_func import read_data_from_file_lines
from thermo.chemical import Chemical
import openpyxl

if __name__ == "__main__":

    # Соотношение R-22 к R-115
    # Общее соотношения газа к маслу с учетом распределительного коэффициента
    g = 0.7
    zs = [0.488 * g, 0.512 * g, 1 - g]  # R-22, R-115, ISO 32
    T = read_data_from_file_lines('temperature.txt')
    P = read_data_from_file_lines('pressure.txt')
    fractions = [read_data_from_file_lines('fraction_r22.txt'), read_data_from_file_lines('fraction_r115.txt')]
    Fraction_R22, Fraction_R115 = fractions
    IDs = ['r-22', 'chloropentafluoroethane', 'water']

    # Молярная масса хладагента R-22
    r22 = Chemical(IDs[0])
    M_r22 = round(r22.MW, 2)  #86.47

    # Молярная масса хладагента R-115
    r115 = Chemical(IDs[1])
    M_r115 = round(r115.MW, 2)  #154.47

    # Критические параметры масла:
    M_oil = 400.0  # Молярная масса масла
    T_oil = 800  # Кельвины
    P_oil = 950e3  # Па
    # Расчет ацентрического фактора по корреляции Ли и Кесслера
    omega = calculate_acentrick_factor(T_oil, P_oil, M_oil)

    # Расчет мольных долей для каждого экспериментального компонента системы
    mole_fractions_r22, mole_fractions_r115, mole_fractions_oil = mass_fractions_to_mole_fractions(
        Fraction_R22, Fraction_R115, M_r22, M_r115, M_oil
    )

    # Расчет мольных долей для каждого компонента системы реализуемой программы
    zs_r22_mole, zs_r115_mole, zs_oil_mole = mass_fractions_to_mole_fractions(
        [[zs[0]]], [[zs[1]]], M_r22, M_r115, M_oil
    )
    zs_mole = [zs_r22_mole[0][0], zs_r115_mole[0][0], zs_oil_mole[0][0]]

    oil_PRMIX = Oil(M_oil, T_oil, P_oil, omega)
    chemical_constants_PR = create_constants_package(oil_PRMIX, IDs)

    # ================= Оптимизация коэффициентов =================
    # Оптимизируемые коэффициенты для R-22 и R-115
    initial_guess_k13_k23 = [0.21, 0.001, 0.35, 0.001, 0.1, 0.001]
    flash_calculator = FlashCalculator(PRMIX, oil_PRMIX, chemical_constants_PR, initial_guess_k13_k23)
    result = flash_calculator.optimize_parameters(initial_guess_k13_k23, [mole_fractions_r22, mole_fractions_r115],
                                                  T, P, zs_mole)

    # Счетчик точек и суммы концентраций
    count = 0
    sum_R22 = 0
    sum_R115 = 0

    workbook = openpyxl.Workbook()
    sheet = workbook.active

    # Записываем заголовки в файл Excel
    sheet.append(
        ['Т(°C)', 'P(кПа)', 'x_R-22(модел.)', 'х_R-115(модел.)', 'х_ISO 32(модел.)', 'x_R-22(эксп.)', 'x_R-115(эксп.)',
         'x_ISO 32(эксп.)',
         'Aбс.ошибка R-22', 'Абс.ошибка R-115'])

    # Цикл расчета вспышки и ошибок
    for i in range(len(T)):
        countT = 0
        sum_R22T = 0
        sum_R115T = 0
        # Массивы данных концентраций газов, давления
        R22_plot = []
        R115_plot = []
        P_bar_plot = []
        Fraction_R22_plot = []
        Fraction_R115_plot = []
        for j in range(len(P[i])):
            # Вызов функции для получения мольных долей, температуры и давления
            flash_data = flash_calculator.calculate_flash(T, P, zs_mole, i, j)
            if flash_data[3] == 1:
                print('================ПРИСУТСТВУЕТ ТОЛЬКО ПАР================')
                print('Температура расчета:', T[i], '°C\n')
                print('Давление расчета:', P[i][j], 'Бар\n')
                print('Соотношение компонентов смеси в газе R-22, R-115, ISO 32:\n', flash_data[0], '\n')
                print('Пар:', flash_data[3], '\n')
                # Обнуляем наши значения долей, если присутствует только пар
                flash_data[0][0] = 0.0
                flash_data[0][1] = 0.0
            else:
                print('================================================================================')
                print('Температура расчета:', T[i], '°C\n')
                print('Давление расчета:', P[i][j], 'Бар\n')
                print('Соотношение компонентов смеси R-22, R-115, ISO 32:\n', flash_data[0], '\n')
                print('Пар:', flash_data[3], '\n')

            # Добавление точек в листы
            R22_plot.append(flash_data[0][0])
            R115_plot.append(flash_data[0][1])
            P_bar_plot.append(P[i][j])
            # Экспериментальные (исходные) данные
            Fraction_R22_plot.append(mole_fractions_r22[i][j])
            Fraction_R115_plot.append(mole_fractions_r115[i][j])

            # Абсолютные ошибки
            ABS_error_R22 = abs(mole_fractions_r22[i][j] - flash_data[0][0])
            ABS_error_R115 = abs(mole_fractions_r115[i][j] - flash_data[0][1])

            # Относительные ошибки
            Rel_error_R22 = ABS_error_R22 / mole_fractions_r22[i][j]
            Rel_error_R115 = ABS_error_R115 / mole_fractions_r115[i][j]

            # Счетчик точек
            count += 1
            countT += 1
            sum_R22T += ABS_error_R22
            sum_R115T += ABS_error_R115

            sum_R22 += ABS_error_R22
            sum_R115 += ABS_error_R115

            print('****** Абсолютная погрешность ******\n')
            print('Погрешность для R-22:', ABS_error_R22, '\n')
            print('Погрешность для R-115:', ABS_error_R115, '\n')
            print('****** Относительная погрешность ******\n')
            print('Погрешность для R-22:', Rel_error_R22, '\n')
            print('Погрешность для R-115:', Rel_error_R115, '\n')

            # Добавляем данные в файл Excel
            sheet.append([
                T[i],
                P[i][j],
                float(flash_data[0][0]),
                float(flash_data[0][1]),
                float(flash_data[0][2]),
                mole_fractions_r22[i][j], mole_fractions_r115[i][j], mole_fractions_oil[i][j],
                float(ABS_error_R22), float(ABS_error_R115)
            ])

        print('****** Значения kijs ******')
        print('k12=', flash_data[4], '\n')
        print('k13=', flash_data[5], '\n')
        print('k23=', flash_data[6], '\n')

        # Средняя ошибка для каждой температуры
        Avg_error_R22T = sum_R22T / countT
        Avg_error_R115T = sum_R115T / countT

        # Построение графика
        plt.title('График зависимости x от P при Т=' + str(T[i]) + '°C', fontsize=18,
                  fontname='Times New Roman')
        plt.xlabel('P, бар', fontsize=18)
        plt.ylabel('x, мольн. доля', fontsize=18)
        plt.scatter(P_bar_plot, R22_plot)
        plt.scatter(P_bar_plot, R115_plot)
        plt.scatter(P_bar_plot, Fraction_R22_plot)
        plt.scatter(P_bar_plot, Fraction_R115_plot)
        plt.plot(P_bar_plot, R22_plot)
        plt.plot(P_bar_plot, R115_plot)
        plt.legend(['R-22 модел.', 'R-115 модел.', 'R-22 эксп.', 'R-115 эксп.'])
        plt.show()

    # Средняя ошибка общая
    Avg_error_R22 = sum_R22 / count
    Avg_error_R115 = sum_R115 / count

    print('****** Средняя погрешность ******')
    print('Погрешность для R-22:', Avg_error_R22, '; В процентах:', round(Avg_error_R22 * 100, 2), '%', '\n')
    print('Погрешность для R-115:', Avg_error_R115, '; В процентах:', round(Avg_error_R115 * 100, 2), '%', '\n')

# Сохраняем файл Excel
workbook.save("result_PR.xlsx")
