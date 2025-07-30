# Пример словаря M_gas_with_oil
M_gas_with_oil = {
    'R-22': 43,
    'R-115': 22,
    'ISO 32': 3
}

# Пример значения flash_data
flash_data = [0.4, 0.3, 0.3]

# Формирование строки соотношения компонентов
components = ', '.join([f'{key}: {value}' for key, value in zip(M_gas_with_oil.keys(), flash_data)])

# Вывод строки
print(f'Соотношение компонентов смеси {components}')