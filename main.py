from matrix import Matrix
from models.specie import Specie
from spatial_ssa import SpatialSSA

matrix: Matrix = Matrix((50, 50, 1))
matrix.put_specie_in_cell(25, 26, 0, 4)
matrix.put_specie_in_cell(0, 0, 0, 12)

ssa: SpatialSSA = SpatialSSA(
    matrix,
    [Specie("H", 0.5, 0)],
    []
)

ssa.draw_animate_plot()
