from random import random, randint

from numpy import shape

from matrix import Matrix
from models.reaction import Reaction
from models.specie import Specie
from spatial_ssa import SpatialSSA

species: list[Specie] = [Specie("A", 1, 0), Specie("B", 1, 1), Specie('C', 1, 2)]
matrix: Matrix = Matrix((50, 50, len(species)))

#
# matrix.put_specie_in_cell(2, 2, 2, 1)

matrix.put_specie_in_cell(25, 22, 0, 10000)
matrix.put_specie_in_cell(25, 27, 1, 10000)



ssa: SpatialSSA = SpatialSSA(
    matrix,
    species,
    [
        Reaction(
            [1, 1, 0],
            [0, 0, 1],
            10,
        ),
        Reaction(
            [1, 0, 0],
            [0, 0, 0],
            1,
        ),
        Reaction(
            [0, 1, 0],
            [0, 0, 0],
            1,
        )
    ]
)

ssa.draw_animate_plot(size=35000, interval=1)
