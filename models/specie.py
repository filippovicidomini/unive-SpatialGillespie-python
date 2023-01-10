class Specie:
    name: str
    diffusion_rate: float
    id: int

    def __init__(self, name: str, diffusion_rate: float, id: int):
        self.name = name
        self.diffusion_rate = diffusion_rate
        self.id = id

    def __str__(self):
        return f"Specie {self.name} with diffusion rate {self.diffusion_rate}"

    def __repr__(self):
        return str(self)
