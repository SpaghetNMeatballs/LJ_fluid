def read_data_from_file_lines(file_path):  # Функция считывания данных из файла
    data = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        if len(lines) == 1:
            data = [float(x) for x in lines[0].strip().split(',')]
        else:
            for line in lines:
                values = line.strip().split(',')
                data.append([float(v) for v in values])
    return data