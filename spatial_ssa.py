from math import floor
from random import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from numpy import ndarray, zeros, shape, log

from event_queue import EventQueue
from matrix import Matrix
from models.reaction import Reaction
from models.specie import Specie
from models.subvolume import SubVolume


class SpatialSSA:
    matrix: Matrix
    subvolume_size: float
    reactions: list[Reaction]
    species: list[Specie]

    subvolumes_diffusion_rates: ndarray
    subvolumes_reaction_rates: ndarray
    subvolumes_next_event_times: ndarray

    def __init__(self, matrix: Matrix, species: list[Specie], reactions: list[Reaction],
                 subvolume_size: float = 1.0):
        self.matrix = matrix
        self.reactions = reactions
        self.species = species
        self.subvolume_size = subvolume_size

        # Data coherence check
        assert (len(self.species) == self.matrix.underlying_matrix.shape[2])

        for reaction in reactions:
            assert (len(reaction.reactants) == len(reaction.products) == len(self.species))

        self.diffusion_rates_sum = sum([specie.diffusion_rate for specie in self.species])

    def recalculate_diffusion_rates_for_subvolume(self, x: int, y: int):
        """
        Recalculates the diffusion rates for a given subvolume.
        :param x: int
        :param y: int
        :return: None
        """
        for specie in self.species:
            self.subvolumes_diffusion_rates[x, y, specie.id] = self.matrix.get_specie_concentration_in_cell(
                x, y,
                specie.id) * specie.diffusion_rate / self.subvolume_size ** 2

    def get_subvolumes_diffusion_rates(self) -> ndarray:
        """
        Returns a 3D array of diffusion rates for each subvolume and specie.
        :return: ndarray
        """

        if self.subvolumes_diffusion_rates is None:
            self.subvolumes_diffusion_rates: ndarray = zeros(self.matrix.underlying_matrix.shape, dtype=float)

            for x in range(shape(self.matrix.underlying_matrix)[0]):
                for y in range(shape(self.matrix.underlying_matrix)[1]):
                    self.recalculate_diffusion_rates_for_subvolume(x, y)

        return self.subvolumes_diffusion_rates

    def get_subvolumes_diffusion_rates_sum(self) -> ndarray:
        """
        Returns a 2D array of diffusion rates for each subvolume.
        :return: ndarray
        """
        return self.get_subvolumes_diffusion_rates().sum(axis=2)

    def _get_h(self, index):
        """
            Calculate the combinations of reactants.
        """
        vector = self._reactants_stoic[index]
        if sum(vector) == 0:
            return 1  # all zeros
        elif sum(vector) == 1:  # just one reactant
            return self._state[np.where(vector == 1)][0]

        # Covering the general case
        else:
            # Case 1
            if sum(vector) == 2 and list(vector).count(1) == 2:
                return np.prod(np.array([self._state[i] for i in np.where(vector == 1)]))

            # Case 2
            if sum(vector) == 2 and 2 in vector:
                return self._state[np.where(vector == 1)] * (self._state[np.where(vector == 1)] - 1) / 2
            # Case 3
            if sum(vector) == 3 and list(vector).count(1) == 3:
                return np.prod(np.array([self._state[i] for i in np.where(vector == 1)]))

            # Case 4
            if sum(vector) == 3 and 1 in vector and 2 in vector:
                return self._state[np.where(vector == 1)] * self._state[np.where(vector == 2)] * (
                            self._state[np.where(vector == 2)] - 1) / 2

            # Case 5
            if sum(vector) == 3 and 3 in vector:
                return self._state[np.where(vector == 3)] * (self._state[np.where(vector == 3)] - 1) * (
                            self._state[np.where(vector == 3)] - 2) / 6



    def get_subvolume_reaction_rate(self, x: int, y: int):
        species_vector = self.matrix[x, y, :] # specie nel sottovolume interessato
        # ora dobbiamo calcolare le propensity per ogni reazione





    def get_subvolumes_reaction_rates(self) -> ndarray:
        """
        Returns a 3D array of reaction rates for each subvolume and specie.
        :return: ndarray
        """

        if self.subvolumes_reaction_rates is None:
            self.subvolumes_reaction_rates: ndarray = zeros(
                (self.matrix.underlying_matrix.shape[0], self.matrix.underlying_matrix.shape[1], len(self.reactions)),
                dtype=float)

            for x in range(shape(self.matrix.underlying_matrix)[0]):
                for y in range(shape(self.matrix.underlying_matrix)[1]):
                    for reaction_id in range(len(self.reactions)):
                        # TODO: Reimplementare!
                        self.subvolumes_reaction_rates[x, y, reaction_id] = 0.0  # self.reactions[reaction_id].rate

        return self.subvolumes_reaction_rates

    def get_subvolumes_reaction_rate_sum(self) -> ndarray:
        """
        Returns a 2D array of reaction rates for each subvolume.
        :return: ndarray
        """
        return self.get_subvolumes_reaction_rates().sum(axis=2)

    def get_subvolumes_total_rate(self) -> ndarray:
        """
        Returns a 2D array of total rates for each subvolume.
        :return: ndarray
        """
        return self.get_subvolumes_diffusion_rates_sum() + self.get_subvolumes_reaction_rate_sum()

    def get_subvolumes_next_event_times(self) -> ndarray:
        """
        Returns the time of the next event.
        :return: float
        """
        if self.subvolumes_next_event_times is None:
            subvolumes_total_rate: ndarray = self.get_subvolumes_total_rate()
            self.subvolumes_next_event_times: ndarray = zeros(subvolumes_total_rate.shape, dtype=float)

            for x in range(shape(subvolumes_total_rate)[0]):
                for y in range(shape(subvolumes_total_rate)[1]):
                    if subvolumes_total_rate[x, y] == 0.0:
                        self.subvolumes_next_event_times[x, y] = None
                    else:
                        self.subvolumes_next_event_times[x, y] = -log(random()) / subvolumes_total_rate[x, y]

        return self.subvolumes_next_event_times

    def clear_cache(self):
        self.subvolumes_diffusion_rates = None
        self.subvolumes_reaction_rates = None
        self.subvolumes_next_event_times = None

    def initialize(self):
        self.clear_cache()

        self.event_queue: EventQueue = EventQueue()

        for x in range(self.get_subvolumes_next_event_times().shape[0]):
            for y in range(self.get_subvolumes_next_event_times().shape[1]):
                if self.get_subvolumes_next_event_times()[x, y] > 0.0:
                    self.event_queue.insert(
                        SubVolume(
                            (x, y),
                            self.get_subvolumes_diffusion_rates()[x, y],
                            self.get_subvolumes_reaction_rates()[x, y],
                            self.get_subvolumes_next_event_times()[x, y]
                        )
                    )

    def step(self, i: int):
        im.set_data(self.matrix.underlying_matrix[:, :, 0])

        next_event: SubVolume = self.event_queue.extract_min()
        rand: float = random()

        if rand < self.get_subvolumes_reaction_rates()[next_event.coordinates] / \
                self.get_subvolumes_total_rate()[next_event.coordinates]:
            # Chemical reaction
            pass
        else:
            # Diffusion
            diffusion_rate_rand: float = rand * self.diffusion_rates_sum
            for specie in sorted(self.species, key=lambda s: s.diffusion_rate):
                if diffusion_rate_rand < specie.diffusion_rate:
                    # Diffusion of specie.id
                    direction: int = floor(rand * 8)
                    moved: bool = False

                    while not moved:
                        try:
                            self.matrix.move_specie(*next_event.coordinates, specie.id, direction)
                            moved = True
                        except IndexError:
                            direction = (direction + 1) % 8

                    self.recalculate_diffusion_rates_for_subvolume(*next_event.coordinates)
                    self.recalculate_diffusion_rates_for_subvolume(*self.matrix.get_neighbour_coordinates(
                        *next_event.coordinates, direction))

                    self.event_queue.remove_with_coordinates(next_event.coordinates)
                    self.event_queue.remove_with_coordinates(self.matrix.get_neighbour_coordinates(
                        *next_event.coordinates, direction))

                    if self.get_subvolumes_total_rate()[next_event.coordinates] > 0.0:
                        self.event_queue.insert(
                            SubVolume(
                                next_event.coordinates,
                                self.get_subvolumes_diffusion_rates()[next_event.coordinates],
                                self.get_subvolumes_reaction_rates()[next_event.coordinates],
                                -log(random()) / self.get_subvolumes_total_rate()[next_event.coordinates],
                            ))

                    if self.get_subvolumes_total_rate()[self.matrix.get_neighbour_coordinates(
                            *next_event.coordinates, direction)]:
                        self.event_queue.insert(
                            SubVolume(
                                self.matrix.get_neighbour_coordinates(*next_event.coordinates, direction),
                                self.get_subvolumes_diffusion_rates()[self.matrix.get_neighbour_coordinates(
                                    *next_event.coordinates, direction)],
                                self.get_subvolumes_reaction_rates()[self.matrix.get_neighbour_coordinates(
                                    *next_event.coordinates, direction)],
                                -log(random()) / self.get_subvolumes_total_rate()[self.matrix.get_neighbour_coordinates(
                                    *next_event.coordinates, direction)],
                            )
                        )

                    break

    def draw_animate_plot(self):
        # Animate a plot with the simulation using matplotlib
        global im, target_loop_count
        target_loop_count = 1000

        fig, ax = plt.subplots()
        im = ax.imshow(self.matrix.underlying_matrix[:, :, 0], cmap="tab20b", interpolation='nearest')

        self.initialize()
        ani = animation.FuncAnimation(fig, self.step, save_count=target_loop_count)
        ani.save('basic_animation.mp4', fps=30)
