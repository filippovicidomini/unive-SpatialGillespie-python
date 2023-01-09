import matplotlib.pyplot as plt
from numpy import ndarray, zeros


class Matrix:
    underlying_matrix: ndarray

    def __init__(self, shape: tuple = (10, 10, 1)):
        self.underlying_matrix = zeros(shape, dtype=float)

    def get_adjacent(self, x: int, y: int) -> ndarray:
        return self.underlying_matrix[x - 1:x + 2, y - 1:y + 2]

    def put_specie_in_cell(self, x: int, y: int, specie_id: int, specie_concentration: float):
        self.underlying_matrix[x, y, specie_id] = specie_concentration

    def get_specie_concentration_in_cell(self, x: int, y: int, specie_id: int) -> float:
        return self.underlying_matrix[x, y, specie_id]

    def plot_concentrations(self, specie_id: int):
        plt.matshow(self.underlying_matrix[:, :, specie_id])
        plt.show()

    def move_specie(self, x: int, y: int, specie_id: int, direction: int):
        if direction == 0:
            self.underlying_matrix[x + 1, y, specie_id] += 1
        elif direction == 1:
            self.underlying_matrix[x, y + 1, specie_id] += 1
        elif direction == 2:
            self.underlying_matrix[x - 1, y, specie_id] += 1
        elif direction == 3:
            self.underlying_matrix[x, y - 1, specie_id] += 1
        elif direction == 4:
            self.underlying_matrix[x + 1, y + 1, specie_id] += 1
        elif direction == 5:
            self.underlying_matrix[x - 1, y + 1, specie_id] += 1
        elif direction == 6:
            self.underlying_matrix[x - 1, y - 1, specie_id] += 1
        elif direction == 7:
            self.underlying_matrix[x + 1, y - 1, specie_id] += 1

        self.underlying_matrix[x, y, specie_id] -= 1

    @staticmethod
    def get_neighbour_coordinates(x: int, y: int, direction: int) -> tuple:
        if direction == 0:
            return x + 1, y
        elif direction == 1:
            return x, y + 1
        elif direction == 2:
            return x - 1, y
        elif direction == 3:
            return x, y - 1
        elif direction == 4:
            return x + 1, y + 1
        elif direction == 5:
            return x - 1, y + 1
        elif direction == 6:
            return x - 1, y - 1
        elif direction == 7:
            return x + 1, y - 1
        else:
            return x, y
