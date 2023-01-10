from matrix import Matrix
from models.reaction import Reaction
from models.specie import Specie
from spatial_ssa import SpatialSSA

species: list[Specie] = [Specie("A", 0.3, 0), Specie("B", 0.2, 1), Specie("C", 0.5, 2)]
matrix: Matrix = Matrix((10, 10, len(species)))

matrix.put_specie_in_cell(4, 4, 0, 60)
matrix.put_specie_in_cell(6, 6, 1, 60)

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

ssa.draw_animate_plot(500)
