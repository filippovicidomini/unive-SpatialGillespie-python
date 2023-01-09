from matrix import Matrix

honeycomb: Matrix = Matrix((100, 100, 1))

honeycomb.put_specie_in_cell(50, 50, 0, 1.0)

honeycomb.plot_concentrations(0)
