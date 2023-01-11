from matrix import Matrix
from models.reaction import Reaction
from models.specie import Specie
from spatial_ssa import SpatialSSA

species: list[Specie] = [Specie("x", 0.1, 0), Specie("y1", 8, 1), Specie('y2', 10, 2)]
matrix: Matrix = Matrix((15, 15, len(species)))

matrix.put_specie_in_cell(5, 5, 0, 100)
matrix.put_specie_in_cell(7, 7, 1, 1000)
matrix.put_specie_in_cell(6, 6, 2, 1000)
matrix.put_specie_in_cell(10, 7, 0, 100)
matrix.put_specie_in_cell(6, 6, 2, 1000)
matrix.put_specie_in_cell(13, 9, 1, 1000)
matrix.put_specie_in_cell(11, 9, 2, 1000)
matrix.put_specie_in_cell(2, 4, 1, 1000)
matrix.put_specie_in_cell(7, 2, 2, 1000)
matrix.put_specie_in_cell(3, 13, 0, 100)
matrix.put_specie_in_cell(13, 3, 2, 1000)
matrix.put_specie_in_cell(1, 8, 0, 100)
matrix.put_specie_in_cell(4, 8, 0, 100)
matrix.put_specie_in_cell(3, 3, 2, 1000)
matrix.put_specie_in_cell(7, 2, 2, 1000)
matrix.put_specie_in_cell(5, 3, 2, 1000)

# matrix.put_specie_in_cell(2, 2, 2, 1)

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
            10,
        ),
        Reaction(
            [0, 0, 1],
            [0, 0, 0],
            10,
        )
    ]
)

ssa.draw_animate_plot(size=100)
