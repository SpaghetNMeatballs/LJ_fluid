from Flashers_and_classes.classes import temperature, pressure, Fraction_R115, Fraction_R22
from thermo import CEOSGas, CEOSLiquid, FlashVL

# Соотношение R-22 к R-115
# Общее соотношения газа к маслу с учетом распределительного коэффициента
g = 0.8
zs = [0.488 * g, 0.512 * g, 1 - g]


# Функция оптимизации средней ошибки
def func_opt(parametrs, class_flash, properties, constants):
    # Счетчик точек и суммы концентраций
    count = 0
    sum_R22 = 0
    sum_R115 = 0
    for i in range(len(temperature.values)):
        for j in range(len(pressure.values[i])):

            # kij значения для газов и масла
            k13 = parametrs[0] + (temperature.values[i] - temperature.values[0]) * parametrs[1]
            k23 = parametrs[2] + (temperature.values[i] - temperature.values[0]) * parametrs[3]
            k12 = parametrs[4] + (temperature.values[i] - temperature.values[0]) * parametrs[5]
            kijs = [[0, k12, k13],
                    [k12, 0, k23],
                    [k13, k23, 0]]

            eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
            gas = CEOSGas(class_flash, eos_kwargs=eos_kwargs)
            liquid = CEOSLiquid(class_flash, eos_kwargs=eos_kwargs)
            flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)

            PT = flasher.flash(T=temperature.values[i] + 273.15, P=pressure.values[i][j] * pow(10, 5), zs=zs)
            if PT.VF != 1:  # Не только пар

                # Абсолютные ошибки
                ABS_error_R22 = abs(Fraction_R22.values[i][j] - PT.liquid0.zs[0])
                ABS_error_R115 = abs(Fraction_R115.values[i][j] - PT.liquid0.zs[1])

            else:
                ABS_error_R22 = abs(Fraction_R22[i][j] - 0.0)
                ABS_error_R115 = abs(Fraction_R115[i][j] - 0.0)

            # Счетчик точек
            count += 1
            sum_R22 += ABS_error_R22
            sum_R115 += ABS_error_R115

    # Средняя ошибка
    Avg_error_R22 = sum_R22 / count
    Avg_error_R115 = sum_R115 / count

    return Avg_error_R115 + Avg_error_R22


def calculate_flash(T, P, fractions, i, j, result, class_flash, constants, properties):  # Функция возврата мольных долей, температуры и давления
    k13 = result.x[0] + (T[i] - T[0]) * result.x[1]
    k23 = result.x[2] + (T[i] - T[0]) * result.x[3]
    k12 = result.x[4] + (T[i] - T[0]) * result.x[5]

    kijs = [[0, k12, k13],
            [k12, 0, k23],
            [k13, k23, 0]]

    eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
    gas = CEOSGas(class_flash, eos_kwargs=eos_kwargs)
    liquid = CEOSLiquid(class_flash, eos_kwargs=eos_kwargs)
    flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)

    PT = flasher.flash(T=T[i] + 273.15, P=P[i][j] * pow(10, 5), zs=fractions)

    if PT.VF == 1:  # Есть только пар
        return [PT.gas.zs, T, P, PT.VF, k12, k13, k23]
    else:  # Не только пар
        return [PT.liquid0.zs, T, P, PT.VF, k12, k13, k23]

