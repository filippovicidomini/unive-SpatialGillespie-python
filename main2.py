from random import random, randint

from numpy import shape

from matrix import Matrix
from models.reaction import Reaction
from models.specie import Specie
from spatial_ssa import SpatialSSA

species: list[Specie] = [Specie("A", 1, 0),       # definisco le specie chimiche che interagiscono
                         Specie("B", 1, 1),       # ognuna con nome, costante diffusiva e indice
                         Specie('C', 1, 2)]

matrix: Matrix = Matrix((50, 50, len(species)))   # definisco la matrice di dimensione 50x50x3
                                                  # dove 50x50 sono le dimensioni della griglia
                                                  # e 3 è il numero di specie chimiche che interagiscono


matrix.put_specie_in_cell(25, 25, 0, 100)         # metto 100 molecole di A nella cella (25,25)
matrix.put_specie_in_cell(20, 20, 1, 100)         # metto 100 molecole di B nella cella (20,20)



ssa: SpatialSSA = SpatialSSA(                     # definisco l'oggetto che simula la reazione
    matrix,                                       # passo la matrice
    species,                                      # passo le specie chimiche
    [
        Reaction(                                 # definisco la reazione A+B -> C
            [1, 1, 0],                            # 1 molecola di A e 1 molecola di B
            [0, 0, 1],                            # 1 molecola di C
            10,                                   # definisco la costante di reazione
        )
    ]
)

ssa.draw_animate_plot(size=1000, interval=0.0001) # disegno l'animazione della reazione
                                                  # size è il numero di passi della simulazione
                                                  # interval è il tempo tra un passo e l'altro


