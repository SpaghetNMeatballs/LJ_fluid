from os import listdir


def load_from_directory(dir_path: str, column: int) -> dict[float:dict[float:float]]:
    result = {}
    for file in listdir(dir_path):
        with open(dir_path + "\\" + file, "r") as f:
            text = f.read()
            # Считываем давление из имени
            p = float(file[:-4])
            result[p] = {}
            for line in text.split('\n'):
                cur_line = line.split(' ')
                # Забираем значение из искомой колонки
                t = float(cur_line[0])
                val = float(cur_line[column])
                result[p][t] = val
    return result


class ViscData:
    def __init__(self, dir_path: str):
        self.density = load_from_directory(dir_path, 2)
        self.visc_k = load_from_directory(dir_path, 3)
        self.visc_d = load_from_directory(dir_path, 4)
        self.p = list(self.density.keys())

    def get_by_pressure(self, dict_index: int, p: float) -> list[list[float]]:
        if p not in self.p:
            raise KeyError(f'{p} not in pressure array')
        if not 0 <= dict_index <= 2:
            raise KeyError(f'{dict_index} not in range 0-2')
        val_dict = [self.density, self.visc_k, self.visc_d][dict_index]
        t_val_dict = val_dict[p]
        t_arr = list(t_val_dict.keys())
        return [[p, t, t_val_dict[t]] for t in t_arr]

    def get_all_points(self, dict_index: int) -> list[list[float]]:
        if not 0 <= dict_index <= 2:
            raise KeyError(f'{dict_index} not in range 0-2')

        result = []
        for p in self.p:
            t_val_dict = [self.density, self.visc_k, self.visc_d][dict_index]
            t_arr = list(t_val_dict.keys())
            result += [[p, t, t_val_dict[t]] for t in t_arr]
        return result


def aad(x_exp, x_calc):
    assert (len(x_exp) == len(x_calc))
    returnable = 0
    for i in range(len(x_exp)):
        returnable += abs((x_exp[i] - x_calc[i]) / x_exp[i])
    returnable *= 100 / len(x_exp)
    return returnable


def aad_new(params: list[float], func: callable, points: list, exp: list[float],
            static_params: list[float] = None) -> float:
    if len(points) != len(exp):
        raise ValueError(f"Input lists have different shapes: {len(points)} and {len(exp)}")
    result = 0
    for i in range(len(points)):
        if static_params:
            result += abs((func(params, static_params, points[i]) - exp[i]) / exp[i])
        else:
            result += abs((func(params, points[i]) - exp[i]) / exp[i])
    result *= 100 / len(points)
    return result
