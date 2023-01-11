from matrix import Matrix
from models.reaction import Reaction
from models.specie import Specie
from spatial_ssa import SpatialSSA
from numpy import shape

species: list[Specie] = [Specie("A", 0.5, 0), Specie("B", 0.5, 1), Specie('P', 0.5, 2), Specie('X', 0.5, 3), Specie('Y', 0.5, 4), Specie('Z', 0.5, 5)]
matrix: Matrix = Matrix((50, 50, len(species)))

# EXAMPLE OF HOW PUT MOLECULA IN SUBVOLUME
# matrix.put_specie_in_cell(2, 2, 2, 1)

# put every speces in every subvolume
for i in range(shape(matrix.underlying_matrix)[0]):
    for j in range(shape(matrix.underlying_matrix)[1]):
        for k in range(shape(matrix.underlying_matrix)[2]):
            if k == 2:
                continue

            if k == 4:
                continue
            if k == 3:
                continue
            matrix.put_specie_in_cell(i, j, k, 10)


matrix.put_specie_in_cell(25, 25, 4, 1000)
matrix.put_specie_in_cell(25, 25, 3, 1000)

ssa: SpatialSSA = SpatialSSA(
    matrix,
    species,
    [
        Reaction(
            [1, 0, 0, 0, 1, 0],
            [0, 0, 1, 1, 0, 0],
            2,
        ),
        Reaction(
            [0, 0, 0, 1, 1, 0],
            [0, 0, 2, 0, 0, 0],
            1.8 * (10**6),
        ),
        Reaction(
            [1, 0, 0, 1, 0, 0],
            [0, 0, 0, 2, 0, 2],
            48,
        ),
        Reaction(
            [0, 0, 0, 2, 0, 0],
            [1, 0, 1, 0, 0, 0],
            3 * (10**3),
        ),
        Reaction(
            [0, 2, 0, 0, 0, 2],
            [0, 0, 0, 0, 1, 0],
            1.6,
        )
    ]
)

ssa.draw_animate_plot(size=100)
