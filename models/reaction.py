class Reaction:
    reactants: list[int]
    products: list[int]
    rate: float

    def __init__(self, reactants: list[int], products: list[int], rate: float):
        self.reactants = reactants
        self.products = products
        self.rate = rate

    def __str__(self):
        return f"Reaction(reactants={self.reactants}, products={self.products}, rate={self.rate})"

    def get_h(self, concentrations: list[int]) -> float:
        if sum(self.reactants) == 0:
            return 1
        elif sum(self.reactants) == 1 and 1 in self.reactants:
            return concentrations[self.reactants.index(1)]
        elif sum(self.reactants) == 2 and self.reactants.count(1) == 2:
            support: list = []
            for i in range(len(self.reactants)):
                if self.reactants[i] == 1:
                    support.append(i)
            return concentrations[support[0]] * concentrations[support[1]]
        elif sum(self.reactants) == 2 and 2 in self.reactants:
            return concentrations[self.reactants.index(2)] * (concentrations[self.reactants.index(2)] - 1) / 2
        elif sum(self.reactants) == 3 and self.reactants.count(1) == 3:
            support: list = []
            for i in range(len(self.reactants)):
                if self.reactants[i] == 1:
                    support.append(i)
            return concentrations[support[0]] * concentrations[support[1]] * concentrations[support[2]]
        elif sum(self.reactants) == 3 and 1 in self.reactants and 2 in self.reactants:
            return concentrations[self.reactants.index(1)] * concentrations[self.reactants.index(2)] * (
                    concentrations[self.reactants.index(2)] - 1) / 2
        elif sum(self.reactants) == 3 and 3 in self.reactants:
            return concentrations[self.reactants.index(3)] * (concentrations[self.reactants.index(3)] - 1) * (
                    concentrations[self.reactants.index(3)] - 2) / 6
        elif sum(self.reactants) == 4 and self.reactants.count(2) == 2:
            support: list = []
            for i in range(len(self.reactants)):
                if self.reactants[i] == 2:
                    support.append(i)
            return concentrations[support[0]] * (concentrations[support[0]] - 1) / 2 * concentrations[support[1]] * (
                    concentrations[support[1]] - 1) / 2

    def get_rate(self, concentrations: list[int]) -> float:
        return self.rate * self.get_h(concentrations)