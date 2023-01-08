class SubVolume:
    coordinates: tuple
    diffusion_rate: float
    reaction_rates: list[float]
    next_event_time: float

    def __init__(self, coordinates: tuple, diffusion_rate: float, reaction_rates: list[float], next_event_time: float):
        self.coordinates = coordinates
        self.diffusion_rate = diffusion_rate
        self.reaction_rates = reaction_rates
        self.next_event_time = next_event_time
