from random import random, randint

from numpy import shape

from matrix import Matrix
from models.reaction import Reaction
from models.specie import Specie
from spatial_ssa import SpatialSSA

species: list[Specie] = [Specie("food", 1, 0), Specie("conigli", 10, 1), Specie('wolf', 10, 2)]
matrix: Matrix = Matrix((50, 50, len(species)))

#
# matrix.put_specie_in_cell(2, 2, 2, 1)

# put food in the invirorment
n_food: int = randint(100, 250)
for _ in range(n_food):
    x: int = randint(0, 49)
    y: int = randint(0, 49)
    matrix.put_specie_in_cell(x, y, 0, 10)

# put conigli fam
n_conigli: int = randint(2000, 2500)
for _ in range(n_conigli):
    x: int = randint(0, 49)
    y: int = randint(0, 49)
    matrix.put_specie_in_cell(x, y, 1, 10)

# put wolf fam
n_wolf: int = randint(2000, 2500)
for _ in range(n_wolf):
    x: int = randint(0, 49)
    y: int = randint(0, 49)
    matrix.put_specie_in_cell(x, y, 2, 10)



matrix.put_specie_in_cell(2, 2, 0, 1)

ssa: SpatialSSA = SpatialSSA(
    matrix,
    species,
    [
        Reaction(
            [1, 1, 0],
            [1, 2, 0],
            10,
        ),
        Reaction(
            [0, 1, 1],
            [0, 0, 2],
            0.01,
        ),
        Reaction(
            [0, 0, 1],
            [0, 0, 0],
            10,
        )
    ]
)

ssa.draw_animate_plot(size=10000, interval=1)
