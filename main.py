from matrix import Matrix
from models.reaction import Reaction
from models.specie import Specie
from spatial_ssa import SpatialSSA

species: list[Specie] = [Specie("A", 0.3, 0), Specie("B", 0.2, 1), Specie("C", 0.5, 2)]
matrix: Matrix = Matrix((5, 5, len(species)))

matrix.put_specie_in_cell(0, 0, 0, 60)
matrix.put_specie_in_cell(1, 1, 1, 60)
matrix.put_specie_in_cell(2, 2, 2, 1)

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
