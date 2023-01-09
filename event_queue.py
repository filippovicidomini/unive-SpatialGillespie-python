from models.subvolume import SubVolume


class EventQueue:
    heap: list[SubVolume] = []

    def __init__(self, subvolumes: list[SubVolume] = []):
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

    def heapify_down(self, index: int):
        while index < len(self.heap):
            left_child_index = 2 * index + 1
            right_child_index = 2 * index + 2
            min_child_index = index

            if left_child_index < len(self.heap) and self.heap[left_child_index].next_event_time < self.heap[
                min_child_index].next_event_time:
                min_child_index = left_child_index

            if right_child_index < len(self.heap) and self.heap[right_child_index].next_event_time < self.heap[
                min_child_index].next_event_time:
                min_child_index = right_child_index

            if min_child_index == index:
                break

            self.heap[index], self.heap[min_child_index] = self.heap[min_child_index], self.heap[index]
            index = min_child_index

    def get_min(self) -> SubVolume:
        return self.heap[0]

    def extract_min(self) -> SubVolume:
        min = self.get_min()
        self.heap[0] = self.heap[-1]
        self.heap.pop()
        self.heapify_down(0)
        return min

    def event_count(self) -> int:
        return len(self.heap)

    def remove_element(self, subvolume: SubVolume):
        self.heap.remove(subvolume)
        self.heapify_down(0)

    def remove_with_coordinates(self, coordinates: tuple[int, int]):
        for subvolume in self.heap:
            if subvolume.coordinates == coordinates:
                self.remove_element(subvolume)
                break
