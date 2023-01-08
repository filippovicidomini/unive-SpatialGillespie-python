class Reaction:
    reactants: list[int]
    products: list[int]
    rate: float

    def __init__(self, reactants: list[int], products: list[int], rate: float):
        self.reactants = reactants
        self.products = products
        self.rate = rate
