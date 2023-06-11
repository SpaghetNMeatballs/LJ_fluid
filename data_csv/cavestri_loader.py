import csv
import os
from dataclasses import dataclass


@dataclass
class Point:
    t: float
    p: float
    x: list[float]
    mu: float
    rho: float


class DataManager:
    def __init__(self):
        self.points = []

    def add_point(self, values: list[float]) -> None:
        assert len(values) == 8, 'Wrong len of input list'
        t = values[0]
        p = values[1]
        liq_fractions = values[4:6]
        visc = values[6]
        rho = values[7]
        self.points.append(Point(t, p, liq_fractions, visc, rho))

    def get_by_t(self, t: float) -> list[Point]:
        result = []
        for point in self.points:
            if point.t == t:
                result.append(point)
        return result

    def get_by_p(self, p: float) -> list[Point]:
        result = []
        for point in self.points:
            if point.p == p:
                result.append(point)
        return result


def get_data(part: bool) -> DataManager:
    result = DataManager()
    goal = 0.488 if part else 1
    file_path = os.path.join(os.path.dirname(__file__), 'cavestri_data.csv')
    with open(file_path, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=';')
        for row in reader:
            float_row = [float(i) for i in row]
            if float_row[2] == goal:
                result.add_point(float_row)
    return result


if __name__ == "__main__":
    test = get_data(True)
    print(test.points[0])
