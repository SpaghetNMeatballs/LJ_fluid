import math


# Функция перевода из массовых долей в мольные доли для двух/трехкомпонентной системы с маслом
def mass_fractions_to_mole_fractions(Fraction_1, Fraction_2, M_1, M_2, M_oil):
    mole_fractions_1 = []
    mole_fractions_2 = []
    mole_fractions_oil = []

    for i in range(len(Fraction_1)):
        mole_fractions_row_1 = []
        mole_fractions_row_2 = []
        mole_fractions_row_oil = []

        for j in range(len(Fraction_1[i])):
            fraction_1 = Fraction_1[i][j]

            if Fraction_2:
                fraction_2 = Fraction_2[i][j]
                n_2 = fraction_2 / M_2
            else:
                fraction_2 = 0
                n_2 = 0

            fraction_oil = 1 - (fraction_1 + fraction_2)

            n_1 = fraction_1 / M_1
            n_oil = fraction_oil / M_oil

            total_moles = n_1 + n_2 + n_oil

            mole_fractions_row_1.append(n_1 / total_moles)
            mole_fractions_row_2.append(n_2 / total_moles)
            mole_fractions_row_oil.append(n_oil / total_moles)

        mole_fractions_1.append(mole_fractions_row_1)
        mole_fractions_2.append(mole_fractions_row_2)
        mole_fractions_oil.append(mole_fractions_row_oil)

    if Fraction_2:
        return mole_fractions_1, mole_fractions_2, mole_fractions_oil
    else:
        return mole_fractions_1, mole_fractions_oil


 # Расчет ацентрического фактора по корреляции Ли и Кесслера
def calculate_acentrick_factor(T, P, M):
    T_r = T / 647.3
    P_r = P / 220.6e5

    a = -0.2127 * T_r + 0.2733 * (T_r ** 2) + 1.466 * P_r - 0.4686 * P_r * T_r - 0.691 * (P_r ** 2)
    b = 0.2607 * T_r - 1.278 * (T_r ** 2) + 0.9375 * P_r + 1.01 * P_r * T_r - 0.9479 * (P_r ** 2)
    c = 0.01354 * T_r + 0.06052 * (T_r ** 2) - 0.04536 * P_r - 0.015 * P_r * T_r + 0.3557 * (P_r ** 2)
    omega = math.exp(a + b * M ** (1 / 3) + c * M ** (2 / 3))

    return omega
