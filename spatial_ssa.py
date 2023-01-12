from math import floor
from random import random

import matplotlib.pyplot as plt
import numpy
from matplotlib import animation
from matplotlib.colors import LogNorm
from numpy import ndarray, zeros, shape, log

from event_queue import EventQueue
from matrix import Matrix
from models.reaction import Reaction
from models.specie import Specie
from models.subvolume import SubVolume

cmaps = [
    "Greys",
    "Purples",
    "Blues",
    "Greens",
    "Oranges"
]


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

        self.event_queue: EventQueue = EventQueue()

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

    def update_subvolume_reaction_rates(self, x: int, y: int):
        """
        Updates the reaction rates for a given subvolume.
        :param x: int
        :param y: int
        :return: None
        """
        species_vector: list = list(self.matrix.underlying_matrix[x, y, :])
        reaction_rates: ndarray = numpy.array(
            list(map(lambda reaction: reaction.get_rate(species_vector), self.reactions)))

        self.subvolumes_reaction_rates[x, y] = reaction_rates

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
                    self.update_subvolume_reaction_rates(x, y)

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

    stopped: bool = False

    def step(self, loop_count: int):
        print(f"{loop_count} : {self.matrix.get_total_molecules_count()}")

        for i in range(len(self.species)):
            ims[i].set_data(self.matrix.underlying_matrix[:, :, i])

        next_event: SubVolume = self.event_queue.extract_min()

        if next_event is None:
            if not self.stopped:
                self.stopped = True
            return

        rand: float = random()

        if rand < self.get_subvolumes_reaction_rate_sum()[next_event.coordinates] / self.get_subvolumes_total_rate()[next_event.coordinates]:
            # Reaction
            reaction_rate_rand: float = rand * sum(self.get_subvolumes_reaction_rates()[next_event.coordinates])
            for reaction in sorted(zip(self.reactions, self.get_subvolumes_reaction_rates()[next_event.coordinates]),
                                   key=lambda x: x[1], reverse=True):
                if reaction_rate_rand - reaction[1] < 0:
                    self.matrix.execute_reaction(*next_event.coordinates, reaction[0])
                    break

                reaction_rate_rand -= reaction[1]

            self.recalculate_diffusion_rates_for_subvolume(*next_event.coordinates)
            self.update_subvolume_reaction_rates(*next_event.coordinates)

            if self.get_subvolumes_total_rate()[next_event.coordinates] > 0.0:
                self.event_queue.insert(
                    SubVolume(
                        next_event.coordinates,
                        self.get_subvolumes_diffusion_rates()[next_event.coordinates],
                        self.get_subvolumes_reaction_rates()[next_event.coordinates],
                        -log(random()) / self.get_subvolumes_total_rate()[
                            next_event.coordinates] + next_event.next_event_time
                    )
                )

        else:
            # Diffusion
            diffusion_rate_rand: float = rand * self.get_subvolumes_diffusion_rates_sum()[next_event.coordinates]

            for specie in sorted(self.species,
                                 key=lambda s: s.diffusion_rate * self.matrix.get_specie_concentration_in_cell(
                                     *next_event.coordinates, s.id), reverse=True):
                if diffusion_rate_rand - specie.diffusion_rate * self.matrix.get_specie_concentration_in_cell(
                        *next_event.coordinates, specie.id) < 0:
                    # Diffusion of specie.id
                    direction: int = floor(random() * 8)
                    moved: bool = False

                    while not moved:
                        if self.matrix.can_move(*next_event.coordinates, direction):
                            self.matrix.move_specie(*next_event.coordinates, specie.id, direction)
                            moved = True
                        else:
                            direction = (direction + 1) % 8

                    self.recalculate_diffusion_rates_for_subvolume(*next_event.coordinates)
                    self.recalculate_diffusion_rates_for_subvolume(*self.matrix.get_neighbour_coordinates(
                        *next_event.coordinates, direction))

                    self.update_subvolume_reaction_rates(*next_event.coordinates)
                    self.update_subvolume_reaction_rates(*self.matrix.get_neighbour_coordinates(
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
                                -log(random()) / self.get_subvolumes_total_rate()[
                                    next_event.coordinates] + next_event.next_event_time
                            ))

                    if self.get_subvolumes_total_rate()[self.matrix.get_neighbour_coordinates(
                            *next_event.coordinates, direction)] > 0.0:
                        self.event_queue.insert(
                            SubVolume(
                                self.matrix.get_neighbour_coordinates(*next_event.coordinates, direction),
                                self.get_subvolumes_diffusion_rates()[self.matrix.get_neighbour_coordinates(
                                    *next_event.coordinates, direction)],
                                self.get_subvolumes_reaction_rates()[self.matrix.get_neighbour_coordinates(
                                    *next_event.coordinates, direction)],
                                -log(random()) / self.get_subvolumes_total_rate()[self.matrix.get_neighbour_coordinates(
                                    *next_event.coordinates, direction)] + next_event.next_event_time
                            )
                        )

                    break
                diffusion_rate_rand -= specie.diffusion_rate
        return ims

    def draw_animate_plot(self, size: int = 1000, interval: int = 1):
        # i want a subplot for every species
        global ims, ax

        fig, ax = plt.subplots(1, len(self.species), figsize=(len(self.species) * 3, 3))
        ims = []

        for nd, ax in enumerate(ax.flatten()):
            l = ax.imshow(self.matrix.underlying_matrix[:, :, nd], norm=LogNorm(vmin=1, vmax=1000), cmap='inferno')
            ax.axis('off')
            ax.set_title(self.species[nd].name)
            ims.append(l)

        self.initialize()

        ani = animation.FuncAnimation(fig, self.step, frames=size, interval=interval, blit=False, repeat=False)

        plt.show()

        ani.save('test.mp4', fps=30)

