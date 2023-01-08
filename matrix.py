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
