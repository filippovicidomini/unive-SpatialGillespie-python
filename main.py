from matrix import Matrix
from models.reaction import Reaction
from models.specie import Specie
from spatial_ssa import SpatialSSA

species: list[Specie] = [Specie("X", 0.1, 0), Specie("Y", 1, 1)]
matrix: Matrix = Matrix((60, 60, len(species)))

# EXAMPLE OF HOW PUT MOLECULA IN SUBVOLUME
# matrix.put_specie_in_cell(2, 2, 2, 1)

for i in range(0, 60):
    for j in range(0, 60):
        matrix.put_specie_in_cell(i, j, 0, 10)
matrix.put_specie_in_cell(30, 30, 1, 10000)

ssa: SpatialSSA = SpatialSSA(
    matrix,
    species,
    [
        Reaction(
            [1, 1],
            [0, 2],
            10000,
        ),
        Reaction(
            [0, 1],
            [1000, 0],
            0.1,
        )
    ]
)

ssa.draw_animate_plot(size=1000)
