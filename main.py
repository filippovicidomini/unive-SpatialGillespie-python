from matrix import Matrix
from models.reaction import Reaction
from models.specie import Specie
from spatial_ssa import SpatialSSA

species: list[Specie] = [Specie("A", 0.7, 0), Specie("B", 0.1, 1), Specie("C", 0.7, 2)]
matrix: Matrix = Matrix((15, 15, len(species)))

matrix.put_specie_in_cell(7, 7, 0, 100)
for i in range(0, 15):
    for j in range(0, 15):
        matrix.put_specie_in_cell(i, j, 1, 10)



# matrix.put_specie_in_cell(2, 2, 2, 1)

ssa: SpatialSSA = SpatialSSA(
    matrix,
    species,
    [
        Reaction(
            [1, 1, 0],
            [0, 0, 1],
            0.8,
        )
    ]
)

ssa.draw_animate_plot(size=100)
