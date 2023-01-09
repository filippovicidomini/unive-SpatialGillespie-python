from math import floor
from random import random

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

        # Data coherence check
        assert (len(self.species) == self.matrix.underlying_matrix.shape)

        for reaction in reactions:
            assert (len(reaction.reactants) == len(reaction.products) == len(self.species))

        self.diffusion_rates_sum = sum([specie.diffusion_rate for specie in self.species])

    def get_subvolumes_diffusion_rates(self) -> ndarray:
        """
        Returns a 3D array of diffusion rates for each subvolume and specie.
        :return: ndarray
        """

        if self.subvolumes_diffusion_rates is None:
            self.subvolumes_diffusion_rates: ndarray = zeros(self.matrix.underlying_matrix.shape, dtype=float)

            for x in range(shape(self.matrix.underlying_matrix)[0]):
                for y in range(shape(self.matrix.underlying_matrix)[1]):
                    for specie in self.species:
                        self.subvolumes_diffusion_rates[x, y, specie.id] = self.matrix.get_specie_concentration_in_cell(
                            x, y,
                            specie.id) * specie.diffusion_rate / self.subvolume_size ** 2

        return self.subvolumes_diffusion_rates

    def get_subvolumes_diffusion_rates_sum(self) -> ndarray:
        """
        Returns a 2D array of diffusion rates for each subvolume.
        :return: ndarray
        """
        return self.get_subvolumes_diffusion_rates().sum(axis=2)

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

    def step(self):
        self.clear_cache()

        event_queue: EventQueue = EventQueue()
        for x in range(self.get_subvolumes_next_event_times().shape[0]):
            for y in range(self.get_subvolumes_next_event_times().shape[1]):
                if self.get_subvolumes_next_event_times()[x, y] is None:
                    continue
                else:
                    event_queue.insert(
                        SubVolume(
                            (x, y),
                            self.get_subvolumes_diffusion_rates()[x, y],
                            self.get_subvolumes_reaction_rates()[x, y],
                            self.get_subvolumes_next_event_times()[x, y]
                        )
                    )

        # Let the event happen
        while event_queue.event_count() > 0:
            next_event: SubVolume = event_queue.extract_min()
            rand: float = random()

            if rand < self.get_subvolumes_diffusion_rates_sum()[next_event.coordinates] / \
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
                        try:
                            self.matrix.move_specie(*next_event.coordinates, specie.id, direction)
                        except IndexError:
                            direction = (direction + 4) % 8
                            self.matrix.move_specie(*next_event.coordinates, specie.id, direction)
                        break
