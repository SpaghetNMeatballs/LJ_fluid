# Задача. Предсказывание состава жидкой фазы смеси.
# Входные параметры: Температура C, давление Бар (МПа/10) для каждой точки, вязкость cP, плотность г/мл и состав газа
# R-22/R-115 (48.8/51.2) ------> (0.488/0.512)
# Построить зависимость концентрации от давления

import matplotlib.pyplot as plt
import os
import thermo
from thermo import PRMIX, SRKMIX
from Flashers_and_classes.classes import FlashCalculator, create_constants_package, Oil
from functions.chemical_functions import mass_fractions_to_mole_fractions, calculate_acentrick_factor
from functions.read_func import read_data_from_file_lines
from thermo.chemical import Chemical
import openpyxl
import itertools

if __name__ == "__main__":

    # Блок ввода данных (интерфейс)
    print('************************************************************')
    print('         ПРОГРАММА РАСЧЕТА РАСТВОРИМОСТИ ГАЗОВ')
    print('************************************************************\n')
    print('================Выбор уравнения состояния================')
    flash_class = None
    while not flash_class:
        # Получение выбора пользователя для класса PRMIX
        flash_class_name = input("Введите уравнение состояния для расчета (SRKMIX или PRMIX): ")
        if flash_class_name == "SRKMIX":
            flash_class = SRKMIX
        elif flash_class_name == "PRMIX":
            flash_class = PRMIX
        else:
            print("Некорректные данные. Проверьте ввод.")

    # Ввод газов для системы
    # IDs = ['r-22', 'chloropentafluoroethane', 'water']
    print('================Используемые газы================')
    IDs = []
    print("Введите используемые газы в системе. Последним элементом введите 'water', он будет заменен маслом.\n"
          "Для завершения ввода введите 'end':")
    while True:
        element = input("Введите газ/закончите ввод: ")
        if element == "end":
            break
        # Проверка наличия элемента в библиотеке thermo
        try:
            thermo.chemical.Chemical(element)  # Создаем объект Chemical для проверки наличия вещества
        except ValueError:
            print(f"Ошибка: Вещество '{element}' отсутствует в библиотеке, возможно имеет другое название.")
            continue
        IDs.append(element)

    # Общее соотношения газа к маслу с учетом распределительного коэффициента
    # g = 0.7
    # zs = [0.488 * g, 0.512 * g, 1 - g]  # R-22, R-115, ISO 32
    g = float(input("Введите распределительный коэффициент g: "))
    zs = {}
    for i in range(len(IDs) - 1):
        while True:
            component_ratio = float(input(f"Введите соотношение газа в изначальном состоянии {IDs[i]}: "))
            if 0 <= component_ratio <= 1:
                zs[IDs[i]] = [[component_ratio * g]]
                break
            else:
                print("Некорректное соотношение газа. Пожалуйста, введите значение от 0 до 1.")
    # print(zs)

    # Ввод экспериментальных данных
    print('================Данные эксперимента================')
    temperature_file = input("Введите название файла (путь) с температурой: ")
    while not os.path.exists(temperature_file):
        print("Файл с температурой не найден.")
        temperature_file = input("Введите корректное название файла (путь) с температурой: ")

    T = read_data_from_file_lines(temperature_file)
    t_min = T[0]  # Начальная температура для оптимизации

    pressure_file = input("Введите название файла (путь) с давлением: ")
    while not os.path.exists(pressure_file):
        print("Файл с давлением не найден.")
        pressure_file = input("Введите корректное название файла с давлением: ")

    P = read_data_from_file_lines(pressure_file)
    # print(P)

    fractions = {}
    for i in range(len(IDs) - 1):
        while True:
            file_path = input(f"Введите название файла (путь) с долями газа {IDs[i]}: ")
            if not os.path.exists(file_path):
                print(f"Ошибка: Файл '{file_path}' не существует.")
            else:
                fraction_data = read_data_from_file_lines(file_path)
                fractions[IDs[i]] = fraction_data
                break
    # print(fractions)

    # Критические параметры масла:
    print('================Критические параметры масла================')
    while True:
        try:
            oil_name = str(input("Введите название масла: "))
            break
        except ValueError:
            print("Ошибка: Некорректный ввод. Проверьте правильность вводимого значения.")

    while True:
        try:
            M_oil = float(input("Введите молярную массу масла: "))
            break
        except ValueError:
            print("Ошибка: Некорректный ввод. Проверьте правильность вводимого значения.")

    while True:
        try:
            T_oil = float(input("Введите температуру масла в Кельвинах: "))
            break
        except ValueError:
            print("Ошибка: Некорректный ввод. Проверьте правильность вводимого значения.")

    while True:
        try:
            P_oil = float(input("Введите давление масла в Па: "))
            break
        except ValueError:
            print("Ошибка: Некорректный ввод. Проверьте правильность вводимого значения.")
    '''
    M_oil = 400.0  # Молярная масса масла
    T_oil = 800  # Кельвины
    P_oil = 950e3  # Па
    '''
    # Молярная масса хладагентов
    M_gas = {}
    for i in range(len(IDs) - 1):
        gas = Chemical(IDs[i])
        M_gas[IDs[i]] = round(gas.MW, 2)
    M_gas_with_oil = {**M_gas, 'Oil': M_oil}  # Словарь для газов с маслом

    # Расчет ацентрического фактора по корреляции Ли и Кесслера
    omega = calculate_acentrick_factor(T_oil, P_oil, M_oil)

    # Расчет мольных долей для каждого экспериментального компонента системы
    mole_fractions_experiment = mass_fractions_to_mole_fractions(fractions, M_gas, M_oil)
    # print(mole_fractions_experiment)

    # Расчет мольных долей для каждого компонента системы реализуемой программы
    mole_fractions_model = mass_fractions_to_mole_fractions(zs, M_gas, M_oil)

    # Перевод в одномерный список
    zs_mole = list(itertools.chain.from_iterable(mole_fractions_model.values()))
    zs_mole = flat_list = sum(zs_mole, [])
    # print(zs_mole)

    oil = Oil(M_oil, T_oil, P_oil, omega, oil_name)

    #Перезаписали данные последнего элемента под масло
    chemical_constants = create_constants_package(oil, IDs)

    # ================= Оптимизация коэффициентов =================
    # Оптимизируемые коэффициенты для R-22 и R-115
    initial_guess = [0.21, 0.001, 0.35, 0.001, 0.1, 0.001]
    flash_calculator = FlashCalculator(flash_class, oil, chemical_constants, initial_guess, zs_mole, t_min)
    result = flash_calculator.optimize_parameters(initial_guess, mole_fractions_experiment, T, P, zs_mole)

    workbook = openpyxl.Workbook()
    sheet = workbook.active

    # Создание заголовков (CДЕЛАНО)
    headers = ['Т(°C)', 'P(кПа)'] + [f'x_{element}(модел.)' if element != IDs[-1] else 'x_масло(модел.)' for element in IDs] + \
              [f'x_{element}(эксп.)' if element != IDs[-1] else 'x_масло(эксп.)' for element in IDs] + \
              [f'Aбс.ошибка {element}' for element in IDs[:-1]]
    sheet.append(headers)

    count = 0  # Счетчик точек
    sum_dict = {key: 0 for key in M_gas.keys()}  # Общие суммы ошибок

    # Цикл расчета вспышки и ошибок
    for i in range(len(T)):
        count_current_t = 0
        # Словарь сумм ошибок долей на текущей температуре
        sum_dict_current_t = {key: 0 for key in M_gas.keys()}

        # Средняя ошибка для каждой температуры
        Avg_error_current_t = {key: 0 for key in M_gas.keys()}

        # Массивы данных концентраций газов
        Gases_dict = {}
        for key in M_gas.keys():
            Gases_dict[key] = []

        # Создаем список под давления эксперимента
        P_bar_plot = []
        # Создаем словарь для экспериментальных значений газов эксперимента
        Gases_dict_exp = {}
        for key in M_gas.keys():
            Gases_dict_exp[key] = []

        for j in range(len(P[i])):
            # Вызов функции для получения мольных долей, температуры и давления
            flash_data = flash_calculator.calculate_flash(T[i], P[i][j])
            if flash_data[3] == 1:
                print('================ПРИСУТСТВУЕТ ТОЛЬКО ПАР================')
                print('Температура расчета:', T[i], '°C')
                print('Давление расчета:', P[i][j], 'Бар')
                components = ', '.join([f'{key}: {value}' for key, value in zip(M_gas_with_oil.keys(), flash_data[0])])
                print(f'Соотношение компонентов смеси {components}')
                print('Пар:', flash_data[3], '\n')
                # Обнуляем наши значения долей, если присутствует только пар
                for ind in range(len(flash_data[0])):
                    flash_data[0][ind] = 0
            else:
                print('================================================================================')
                print('Температура расчета:', T[i], '°C')
                print('Давление расчета:', P[i][j], 'Бар')
                components = ', '.join([f'{key}: {value}' for key, value in zip(M_gas_with_oil.keys(), flash_data[0])])
                print(f'Соотношение компонентов смеси {components}')
                print('Пар:', flash_data[3], '\n')

            # Добавление точек модели в словарь
            for idx, key in enumerate(M_gas.keys()):
                Gases_dict[key].append(flash_data[0][idx])

            # Добавляем давление процесса
            P_bar_plot.append(P[i][j])
            # Экспериментальные (исходные) данные
            for key in M_gas.keys():
                Gases_dict_exp[key].append(mole_fractions_experiment[key][i][j])

            # Расчет ошибок
            ABS_error = {key: 0 for key in M_gas.keys()}
            Rel_error = {key: 0 for key in M_gas.keys()}
            for idx, key in enumerate(M_gas.keys()):
                # Абсолютные ошибки
                ABS_error[key] = abs(mole_fractions_experiment[key][i][j] - flash_data[0][idx])
                # Относительные ошибки
                Rel_error[key] = ABS_error[key] / mole_fractions_experiment[key][i][j]

            # Счетчик точек
            count += 1
            count_current_t += 1

            # Расчет сумм абсолютных ошибок по температуре
            for key in M_gas.keys():
                sum_dict_current_t[key] += ABS_error[key]
                sum_dict[key] += sum_dict_current_t[key]

            # Вывод абсолютных ошибок
            print('****** Абсолютная погрешность ******')
            for key in M_gas.keys():
                print(f'Погрешность для {key}: {ABS_error[key]}')

            # Вывод относительных ошибок
            print('****** Относительная погрешность ******')
            for key in M_gas.keys():
                print(f'Погрешность для {key}: {Rel_error[key]}')
            print()

            # Добавляем данные в файл Excel
            sheet.append([
                T[i],
                P[i][j]] +
                [float(flash_data[0][ind]) for ind in range(len(zs_mole))] +
                [mole_fractions_experiment[key][i][j] for key in M_gas_with_oil] +
                [ABS_error[key] for key in M_gas])

        print('****** Cписок значений kijs ******')
        print('****** Значения kijs ******')
        for idx, k in enumerate(flash_data[4], start=1):
            print(f'k{idx} =', k)

        # Средняя ошибка для каждой температуры
        for key in M_gas.keys():
            Avg_error_current_t[key] = sum_dict_current_t[key] / count_current_t

        # Построение графика
        plt.title('График зависимости x от P при Т=' + str(T[i]) + '°C', fontsize=18,
                  fontname='Times New Roman')
        plt.xlabel('P, бар', fontsize=18)
        plt.ylabel('x, мольн. доля', fontsize=18)
        legend_labels = []
        for key in M_gas.keys():
            plt.scatter(P_bar_plot, Gases_dict[key])
            plt.scatter(P_bar_plot, Gases_dict_exp[key])
            plt.plot(P_bar_plot, Gases_dict[key])
            legend_labels.append(f'{key} модел.')
            legend_labels.append(f'{key} эксп.')
        plt.legend(legend_labels)
        plt.show()

    # Средняя ошибка общая
    Avg_error = {key: 0 for key in M_gas.keys()}
    for key in M_gas.keys():
        Avg_error[key] = sum_dict[key] / count

    print('****** Средняя абсолютная погрешность ******')
    for key in M_gas.keys():
        print(f'Погрешность для {key} в процентах:', round(Avg_error[key] * 100, 2), '%')

# Сохранение файла (СДЕЛАНО)
filename = f'{flash_class.__name__}_{len(zs)}_gas.xlsx'.replace("<", "").replace(">", "").replace("'", "")
workbook.save(filename)
