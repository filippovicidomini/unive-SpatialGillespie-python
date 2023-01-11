from matrix import Matrix
from models.reaction import Reaction
from models.specie import Specie
from spatial_ssa import SpatialSSA
from numpy import shape

species: list[Specie] = [Specie("X", 10, 0), Specie("Y", 10, 1)]
matrix: Matrix = Matrix((60, 60, len(species)))

# EXAMPLE OF HOW PUT MOLECULA IN SUBVOLUME
# matrix.put_specie_in_cell(2, 2, 2, 1)

matrix.put_specie_in_cell(20, 20, 0, 1000)
matrix.put_specie_in_cell(40, 40, 1, 1000)

ssa: SpatialSSA = SpatialSSA(
    matrix,
    species,
    [
        Reaction(
            [1, 0],
            [2, 0],
            1,
        ),
        Reaction(
            [1, 1],
            [0, 2],
            10,
        ),
        Reaction(
            [0, 1],
            [0, 0],
            1,
        )

    ]
)

ssa.draw_animate_plot(size=100)
