import os


# Функция считывания данных из файла
def read_data_from_file_lines(file_name: str) -> list[list[float]] | list[float]:
    file_path = os.path.join(os.path.dirname(__file__), '..\\data\\' + file_name)
    with open(file_path, 'r') as f:
        lines = f.readlines()
        data = [[float(v) for v in line.strip().split(',')] for line in lines]
    if len(data) == 1:
        return data[0]
    return data
