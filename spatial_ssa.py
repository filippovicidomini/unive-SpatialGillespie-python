from random import random

from numpy import ndarray, zeros, shape, log

from matrix import Matrix
from models.reaction import Reaction


class SpatialSSA:
    diffusion_rates: ndarray
    matrix: Matrix
    subvolume_size: float
    reactions: list[Reaction]

    def __init__(self, matrix: Matrix, diffusion_rates: ndarray, reactions: list[Reaction],
                 subvolume_size: float = 1.0):
        self.matrix = matrix
        self.diffusion_rates = diffusion_rates
        self.reactions = reactions

        # Data coherence check
        assert (self.diffusion_rates.shape == self.matrix.underlying_matrix.shape)

        for reaction in reactions:
            assert (len(reaction.reactants) == len(reaction.products) == len(self.diffusion_rates.shape))

    def get_subvolumes_diffusion_rate(self) -> ndarray:
        """
        Returns a 3D array of diffusion rates for each subvolume and specie.
        :return: ndarray
        """
        diffusion_rates_matrix: ndarray = zeros(self.matrix.underlying_matrix.shape, dtype=float)

        for x in range(shape(self.matrix.underlying_matrix)[0]):
            for y in range(shape(self.matrix.underlying_matrix)[1]):
                for specie_id in range(shape(self.matrix.underlying_matrix)[2]):
                    diffusion_rates_matrix[x, y, specie_id] = self.matrix.get_specie_concentration_in_cell(x, y,
                                                                                                           specie_id) * \
                                                              self.diffusion_rates[specie_id] / self.subvolume_size ** 2

        return diffusion_rates_matrix

    def get_subvolumes_diffusion_rate_sum(self) -> ndarray:
        """
        Returns a 2D array of diffusion rates for each subvolume.
        :return: ndarray
        """
        return self.get_subvolumes_diffusion_rate().sum(axis=2)

    def get_subvolumes_reaction_rate(self) -> ndarray:
        """
        Returns a 3D array of reaction rates for each subvolume and specie.
        :return: ndarray
        """
        reaction_rate_matrix: ndarray = zeros(
            (self.matrix.underlying_matrix.shape[0], self.matrix.underlying_matrix.shape[1], len(self.reactions)),
            dtype=float)

        for x in range(shape(self.matrix.underlying_matrix)[0]):
            for y in range(shape(self.matrix.underlying_matrix)[1]):
                for reaction_id in range(len(self.reactions)):
                    # TODO: Reimplementare!
                    reaction_rate_matrix[x, y, reaction_id] = 0.0  # self.reactions[reaction_id].rate

        return reaction_rate_matrix

    def get_subvolumes_reaction_rate_sum(self) -> ndarray:
        """
        Returns a 2D array of reaction rates for each subvolume.
        :return: ndarray
        """
        return self.get_subvolumes_reaction_rate().sum(axis=2)

    def get_subvolumes_total_rate(self) -> ndarray:
        """
        Returns a 2D array of total rates for each subvolume.
        :return: ndarray
        """
        return self.get_subvolumes_diffusion_rate_sum() + self.get_subvolumes_reaction_rate_sum()

    def get_subvolumes_next_event_time(self) -> ndarray:
        """
        Returns the time of the next event.
        :return: float
        """
        subvolumes_total_rate: ndarray = self.get_subvolumes_total_rate()
        subvolumes_next_event_time: ndarray = zeros(subvolumes_total_rate.shape, dtype=float)

        for x in range(shape(subvolumes_total_rate)[0]):
            for y in range(shape(subvolumes_total_rate)[1]):
                if subvolumes_total_rate[x, y] == 0.0:
                    subvolumes_next_event_time[x, y] = None

                else:
                    subvolumes_next_event_time[x, y] = -log(random()) / subvolumes_total_rate[x, y]

        return subvolumes_next_event_time
