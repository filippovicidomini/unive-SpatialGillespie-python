from matrix import Matrix
from models.reaction import Reaction
from models.specie import Specie
from spatial_ssa import SpatialSSA


species: list[Specie] = [Specie("A", 0.5, 0), Specie("B", 0.5, 1), Specie("C", 0.5, 2)]
matrix: Matrix = Matrix((15, 15, len(species)))

matrix.put_specie_in_cell(4, 4, 0, 15)
matrix.put_specie_in_cell(6, 6, 1, 15)

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

ssa.draw_animate_plot(1000)
