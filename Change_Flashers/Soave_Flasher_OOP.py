# Задача. Предсказывание состава жидкой фазы смеси.
# Входные параметры: Температура C, давление Бар (МПа/10) для каждой точки, вязкость cP, плотность г/мл и состав газа
# R-22/R-115 (48.8/51.2) ------> (0.488/0.512)
# Построить зависимость концентрации от давления

import matplotlib.pyplot as plt
from thermo import SRKMIX
from Flashers_and_classes.classes import FlashCalculator, create_constants_package, Oil
from functions.read_func import read_data_from_file_lines

if __name__ == "__main__":

    # Соотношение R-22 к R-115
    # Общее соотношения газа к маслу с учетом распределительного коэффициента
    g = 0.8
    zs = [0.488 * g, 0.512 * g, 1 - g]
    T = read_data_from_file_lines('temperature.txt')
    P = read_data_from_file_lines('pressure.txt')
    fractions = [read_data_from_file_lines('fraction_r22.txt'), read_data_from_file_lines('fraction_r115.txt')]
    Fraction_R22, Fraction_R115 = fractions
    IDs = ['r-22', 'chloropentafluoroethane', 'water']
    oil_SRKMIX = Oil(400.0, 800, 950e3, 0.2)
    chemical_constants_SR = create_constants_package(oil_SRKMIX, IDs)

    # ================= Оптимизация коэффициентов =================
    # Оптимизируемые коэффициенты для R-22 и R-115
    initial_guess_k13_k23 = [0.21, 0.001, 0.35, 0.001, 0.1, 0.001]
    flash_calculator = FlashCalculator(SRKMIX, oil_SRKMIX, chemical_constants_SR, initial_guess_k13_k23)
    result = flash_calculator.optimize_parameters(initial_guess_k13_k23, fractions, T, P, zs)

    # Счетчик точек и суммы концентраций
    count = 0
    sum_R22 = 0
    sum_R115 = 0
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
            flash_data = flash_calculator.calculate_flash(T, P, zs, i, j)
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
                print('Температура расчета:', T[i], '°C\n')
                print('Давление расчета:', P[i][j], 'Бар\n')
                print('Соотношение компонентов смеси R-22, R-115, ISO 32:\n', flash_data[0], '\n')
                print('Пар:', flash_data[3], '\n')

            # Добавление точек в листы
            R22_plot.append(flash_data[0][0])
            R115_plot.append(flash_data[0][1])
            P_bar_plot.append(P[i][j])
            # Экспериментальные (исходные) данные
            Fraction_R22_plot.append(Fraction_R22[i][j])
            Fraction_R115_plot.append(Fraction_R115[i][j])

            # Абсолютные ошибки
            ABS_error_R22 = abs(Fraction_R22[i][j] - flash_data[0][0])
            ABS_error_R115 = abs(Fraction_R115[i][j] - flash_data[0][1])

            Rel_error_R22 = ABS_error_R22 / Fraction_R22[i][j]
            Rel_error_R115 = ABS_error_R115 / Fraction_R115[i][j]

            # Счетчик точек
            count += 1
            countT += 1
            sum_R22T += ABS_error_R22
            sum_R115T += ABS_error_R115

            sum_R22 += ABS_error_R22
            sum_R115 += ABS_error_R115

            print('Температура расчета:', T[i], '°C\n')
            print('Давление расчета:', P[i][j], 'Бар\n')
            print('Соотношение компонентов смеси R-22, R-115, ISO 32:\n', flash_data[0], '\n')
            print('Пар:', flash_data[3], '\n')
            print('****** Абсолютная погрешность ******\n')
            print('Погрешность для R-22:', ABS_error_R22, '\n')
            print('Погрешность для R-115:', ABS_error_R115, '\n')
            print('****** Относительная погрешность ******\n')
            print('Погрешность для R-22:', Rel_error_R22, '\n')
            print('Погрешность для R-115:', Rel_error_R115, '\n')

        print('****** Значения kijs ******')
        print('k12=', flash_data[4], '\n')
        print('k13=', flash_data[5], '\n')
        print('k23=', flash_data[6], '\n')

        # Средняя ошибка для каждой температуры
        Avg_error_R22T = sum_R22T / countT
        Avg_error_R115T = sum_R115T / countT

        # Построение графика
        plt.title('График зависимости ω от P при Т=' + str(T[i]) + '°C', fontsize=18,
                  fontname='Times New Roman')
        plt.xlabel('P, кПа', fontsize=18)
        plt.ylabel('ω, масс. доля', fontsize=18)
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
    print('Погрешность для R-22:', Avg_error_R22, '; В процентах:', round(Avg_error_R22 * 100, 2), '%'', \n')
    print('Погрешность для R-115:', Avg_error_R115, '; В процентах:', round(Avg_error_R115 * 100, 2), '%', '\n')