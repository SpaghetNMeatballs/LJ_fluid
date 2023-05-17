import csv


class DataManager:
    def __init__(self):
        self.T = []
        self.P = []
        self.gas_fractions = []
        self.liq_fractions = []
        self.visc = []
        self.rho = []

    def add_point(self, values: list[float]) -> None:
        assert len(values) == 8, 'Wrong len of input list'
        self.T.append(values[0])
        self.P.append(values[1])
        self.gas_fractions.append(values[2:4])
        self.liq_fractions.append(values[4:6])
        self.visc.append(values[6])
        self.rho.append(values[7])


def get_data() -> DataManager:
    result = DataManager()
    with open('cavestri_data.csv', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=';')
        for row in reader:
            float_row = [float(i) for i in row]
            result.add_point(float_row)
    return result


if __name__ == "__main__":
    test = get_data()
    print(test.T)
