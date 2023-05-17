# Задача. Предсказывание состава жидкой фазы смеси.
# Входные параметры: Температура C, давление Бар (МПа/10) для каждой точки, вязкость cP, плотность г/мл и состав газа
# R-22/R-115 (48.8/51.2) ------> (0.488/0.512)
# Построить зависимость концентрации от давления


from vdW2_mix import SRKMIX_vdW2_mixing_rule
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from scipy.optimize import minimize
from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, FlashVL
from thermo.interaction_parameters import IPDB

# R-22, R-115, Oil (Iso32)
constants, properties = ChemicalConstantsPackage.from_IDs(['r-22', 'chloropentafluoroethane', 'water'])

# Критические параметры масла:
M_oil = 400.0  # Молярная масса масла
T_oil = 800  # Кельвины
P_oil = 950e3  # Па
Omega_oil = 0.2  # Константа омега

# Перезаписали критические параметры для воды под масло ISO 32
constants2 = ChemicalConstantsPackage(MWs=constants.MWs[:-1] + [M_oil], names=constants.names[:-1] + ['ISO 32'],
                                      omegas=constants.omegas[:-1] + [Omega_oil], Pcs=constants.Pcs[:-1] + [P_oil],
                                      Tcs=constants.Tcs[:-1] + [T_oil])

# Параметры: массив температур
T = [-10, 0, 20, 40, 70, 80, 125]

# Давление R-22: от -10 градусов до 125
P_bar = [[3.97, 3.8, 3.38, 3.05, 2.81, 2.24],  # -10
         [5.55, 4.85, 4.12, 3.17, 2.28],  # 0
         [10.14, 9.66, 5.14, 2.59],  # 20
         [16.767, 14.352, 10.695, 4.8645],  # 40
         [31.395, 25.3575, 18.699, 10.557, 2.553],  # 70
         [31.119, 25.9095, 7.59, 1.725],  # 80
         [32.20575, 18.94045, 8.142, 2.39755]]  # 125

# Вязкость R-22: от -10 градусов до 125
Viscosity = [[20.286, 26.584, 48.671, 83.596, 123.132, 86.521],  # -10
             [12.905, 24.874, 55.281, 84.659, 101.942],  # 0
             [2.536, 4.706, 20.602, 48.295],  # 20
             [1.613, 36.5, 6.348, 11.599],  # 40
             [1.284, 1.959, 2.927, 4.53, 7.08],  # 70
             [1.863, 2.232, 3.46, 4.263],  # 80
             [1.492, 1.593, 1.877, 1.938]]  # 125

# Плотность R-22:
Density = [[0.9708, 0.9687, 0.9494, 0.9384, 0.9265, 0.9255],  # -10
           [0.9728, 0.9539, 0.9258, 0.9169, 0.8962],  # 0
           [1.0286, 0.9827, 0.9124, 0.8889],  # 20
           [1.0027, 0.9487, 0.9183, 0.8747],  # 40
           [0.9835, 0.9391, 0.9087, 0.8653, 0.8414],  # 70
           [0.9176, 0.8925, 0.8399, 0.8259],  # 80
           [0.8785, 0.8488, 0.8273, 0.8073]]  # 125

# Доли R-22 в жидкой фазе
Fraction_R22 = [[0.14240772, 0.13956069, 0.11894281, 0.09661578, 0.06626904, 0.06917368],
                [0.13752263, 0.097602225, 0.08366001, 0.0695826, 0.04349895],
                [0.21429768, 0.1465556, 0.0542013, 0.02756025],
                [0.214346975, 0.14201664, 0.08427174, 0.03830256],
                [0.19733655, 0.1429248, 0.10825416, 0.05892744, 0.01661184],
                [0.10424721, 0.0714132, 0.03954771, 0.04160995],
                [0.08350189, 0.05926864, 0.0315027, 0.01437675]]

# Доли R-115 в жидкой фазе
Fraction_R115 = [[0.06045228, 0.05346931, 0.04211719, 0.04160422, 0.03917096, 0.03790632],
                 [0.09438737, 0.059947775, 0.04484999, 0.0370574, 0.02723105],
                 [0.13246232, 0.0898244, 0.0301587, 0.01326975],
                 [0.146203025, 0.09156336, 0.06461826, 0.02601744],
                 [0.12696345, 0.0803952, 0.05626584, 0.02995256, 0.01222816],
                 [0.05394279, 0.0361368, 0.02234229, 0.02134005],
                 [0.04320811, 0.02999136, 0.0177973, 0.00737325]]


# Общее соотношения газа к маслу (распределительный коэффициент)
g = 0.9
# Соотношение R-22 к R-115 с маслом
zs = [0.488 * g, 0.512 * g, 1 - g]

# Функция оптимизации средней ошибки
def func_opt(parametrs):
    # Счетчик точек и суммы концентраций
    count = 0
    sum_R22 = 0
    sum_R115 = 0
    for i in range(len(T)):
        for j in range(len(P_bar[i])):
            k13 = parametrs[0] + (T[i] - T[0]) * parametrs[1]
            k23 = parametrs[2] + (T[i] - T[0]) * parametrs[3]
            k12 = parametrs[4] + (T[i]-T[0]) * parametrs[5]
            kijs = [[0, k12, k13],
                    [k12, 0, k23],
                    [k13, k23, 0]]
            '''
            l13 = parametrs[6] + (T[i] - T[0]) * parametrs[7]
            l23 = parametrs[8] + (T[i] - T[0]) * parametrs[9]
            l12 = parametrs[10] + (T[i] - T[0]) * parametrs[11]
            '''
            l13 = l23 = l12 = 0.2
            lijs = [[0, l12, l13],
                    [l12, 0, l23],
                    [l13, l23, 0]]

            eos_kwargs = {'Pcs': constants2.Pcs, 'Tcs': constants2.Tcs, 'omegas': constants2.omegas, 'kijs': kijs, 'lijs': lijs}
            gas = CEOSGas(SRKMIX_vdW2_mixing_rule, eos_kwargs=eos_kwargs)
            liquid = CEOSLiquid(SRKMIX_vdW2_mixing_rule, eos_kwargs=eos_kwargs)
            flasher = FlashVL(constants2, properties, liquid=liquid, gas=gas)

            PT = flasher.flash(T=T[i] + 273.15, P=P_bar[i][j] * pow(10, 5), zs=zs)
            if PT.VF == 1:  # Есть только пар
                print(PT.gas.zs)
                print('Неудача')
            else:  # Не только пар

                # Абсолютные ошибки
                ABS_error_R22 = abs(Fraction_R22[i][j] - PT.liquid0.zs[0])
                ABS_error_R115 = abs(Fraction_R115[i][j] - PT.liquid0.zs[1])

                # Относительные ошибки
                Rel_error_R22 = ABS_error_R22 / PT.liquid0.zs[0]
                Rel_error_R115 = ABS_error_R115 / PT.liquid0.zs[1]

                # Счетчик точек
                count += 1
                sum_R22 += ABS_error_R22
                sum_R115 += ABS_error_R115

    # Средняя ошибка
    Avg_error_R22 = sum_R22 / count
    Avg_error_R115 = sum_R115 / count
    return Avg_error_R115 + Avg_error_R22


# ================= Оптимизация коэффициентов =================
# Оптимизируемые коэффициенты для R-22 и R-115
initial_guess_k13_k23 = [0.21, 0.001, 0.35, 0.001, 0.1, 0.001]

result = minimize(func_opt, initial_guess_k13_k23, bounds=((-1, 1), (0, 0.004), (-1, 1), (0, 0.004), (-1, 1), (0, 0.004)))
if result.success:
    print(result.message)
else:
    raise ValueError(result.message)


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
    for j in range(len(P_bar[i])):
        # Коэффициент kij для двух газов R-22 и R-115

        k13 = result.x[0] + (T[i] - T[0]) * result.x[1]
        k23 = result.x[2] + (T[i] - T[0]) * result.x[3]
        k12 = result.x[4] + (T[i] - T[0]) * result.x[5]

        kijs = [[0, k12, k13],
                [k12, 0, k23],
                [k13, k23, 0]]
        l13 = l23 = l12 = -0.4
        lijs = [[0, l12, l13],
                [l12, 0, l23],
                [l13, l23, 0]]
        eos_kwargs = {'Pcs': constants2.Pcs, 'Tcs': constants2.Tcs, 'omegas': constants2.omegas, 'kijs': kijs, 'lijs': lijs}
        gas = CEOSGas(SRKMIX_vdW2_mixing_rule, eos_kwargs=eos_kwargs)
        liquid = CEOSLiquid(SRKMIX_vdW2_mixing_rule, eos_kwargs=eos_kwargs)
        flasher = FlashVL(constants2, properties, liquid=liquid, gas=gas)

        PT = flasher.flash(T=T[i] + 273.15, P=P_bar[i][j] * pow(10, 5), zs=zs)
        if PT.VF == 1:  # Есть только пар
            print(PT.gas.zs)
            print('Неудача')
        else:  # Не только пар

            # Добавление точек в листы
            R22_plot.append(PT.liquid0.zs[0])
            R115_plot.append(PT.liquid0.zs[1])
            P_bar_plot.append(P_bar[i][j])
            Fraction_R22_plot.append(Fraction_R22[i][j])
            Fraction_R115_plot.append(Fraction_R115[i][j])

            # Абсолютные ошибки
            ABS_error_R22 = abs(Fraction_R22[i][j] - PT.liquid0.zs[0])
            ABS_error_R115 = abs(Fraction_R115[i][j] - PT.liquid0.zs[1])

            # Относительные ошибки
            Rel_error_R22 = ABS_error_R22/PT.liquid0.zs[0]
            Rel_error_R115 = ABS_error_R115 / PT.liquid0.zs[1]

            # Счетчик точек
            count += 1
            countT += 1
            sum_R22T += ABS_error_R22
            sum_R115T += ABS_error_R115


            sum_R22 += ABS_error_R22
            sum_R115 += ABS_error_R115

            print('==============================================================================')
            print('Температура расчета:', T[i], '°C; Давление расчета:', P_bar[i][j], 'Бар')
            print('Соотношение компонентов смеси R-22, R-115, ISO 32:', PT.liquid0.zs)
            print('Пар ', PT.VF, '\n')
            print('****** Абсолютная погрешность ******')
            print('Погрешность для R-22:', ABS_error_R22)
            print('Погрешность для R-115:', ABS_error_R115)
            print('****** Относительная погрешность ******')
            print('Погрешность для R-22:', Rel_error_R22)
            print('Погрешность для R-115:', Rel_error_R115)

    # Средняя ошибка для каждой температуры
    Avg_error_R22T = sum_R22T / countT
    Avg_error_R115T = sum_R115T / countT

    # Построение графика
    plt.title('График зависимости C от P при Т='+str(T[i])+'°C', fontsize=18, fontname='Times New Roman')
    plt.xlabel('P', fontsize=18)
    plt.ylabel('C', fontsize=18)
    plt.scatter(P_bar_plot, R22_plot)
    plt.scatter(P_bar_plot, R115_plot)
    plt.scatter(P_bar_plot, Fraction_R22_plot)
    plt.scatter(P_bar_plot, Fraction_R115_plot)
    plt.plot(P_bar_plot, R22_plot)
    plt.plot(P_bar_plot, R115_plot)
    plt.legend(['R-22 модел.', 'R-115 модел.', 'R-22 эксп.', 'R-115 эксп.'])
    plt.show()

# Средняя ошибка общая
Avg_error_R22 = sum_R22/count
Avg_error_R115 = sum_R115/count

print('****** Средняя погрешность ******')
print('Погрешность для R-22:', Avg_error_R22)
print('Погрешность для R-115:', Avg_error_R115)
