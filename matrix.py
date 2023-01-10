import matplotlib.pyplot as plt
from numpy import ndarray, zeros, array

from models.reaction import Reaction


class Matrix:
    underlying_matrix: ndarray

    def __init__(self, shape: tuple = (10, 10, 1)):
        self.underlying_matrix = zeros(shape, dtype=float)

    def get_adjacent(self, x: int, y: int) -> ndarray:
        return self.underlying_matrix[x - 1:x + 2, y - 1:y + 2]

    def put_specie_in_cell(self, x: int, y: int, specie_id: int, specie_concentration: float):
        self.underlying_matrix[x, y, specie_id] = specie_concentration

    def get_specie_concentration_in_cell(self, x: int, y: int, specie_id: int) -> float:
        if self.underlying_matrix[x, y, specie_id] < 0:
            print(f"PANICOOO -) {x, y, specie_id}")
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

        if self.underlying_matrix[x, y, specie_id] < 0:
            print("PANICOOO dopo diff")

    def can_move(self, x: int, y: int, direction: int) -> bool:
        if direction == 0:
            return x < self.underlying_matrix.shape[0] - 1
        elif direction == 1:
            return y < self.underlying_matrix.shape[1] - 1
        elif direction == 2:
            return x > 0
        elif direction == 3:
            return y > 0
        elif direction == 4:
            return x < self.underlying_matrix.shape[0] - 1 and y < self.underlying_matrix.shape[1] - 1
        elif direction == 5:
            return x > 0 and y < self.underlying_matrix.shape[1] - 1
        elif direction == 6:
            return x > 0 and y > 0
        elif direction == 7:
            return x < self.underlying_matrix.shape[0] - 1 and y > 0

    def execute_reaction(self, x: int, y: int, reaction: Reaction):
        print(f"Executing reaction {reaction} in cell ({x}, {y})")
        print(self.underlying_matrix[x, y, :])
        self.underlying_matrix[x, y] += array(reaction.products) - array(reaction.reactants)

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

    def get_total_molecules_count(self):
        print(any(self.underlying_matrix.flatten() > 0))
        return self.underlying_matrix.sum()
