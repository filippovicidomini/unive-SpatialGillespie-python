from models.subvolume import SubVolume


class EventQueue:
    heap: list[SubVolume]

    def __init__(self, subvolumes: list[SubVolume]):
        for subvolume in subvolumes:
            self.insert(subvolume)

    def insert(self, subvolume: SubVolume):
        self.heap.append(subvolume)
        self.heapify_up(len(self.heap) - 1)

    def heapify_up(self, index: int):
        while index > 0:
            parent_index = (index - 1) // 2
            if self.heap[index].next_event_time < self.heap[parent_index].next_event_time:
                self.heap[index], self.heap[parent_index] = self.heap[parent_index], self.heap[index]
                index = parent_index
            else:
                break
