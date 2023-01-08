import unittest
from random import shuffle

from event_queue import EventQueue
from models.subvolume import SubVolume


class EventQueueTest(unittest.TestCase):
    def test_event_queue(self):
        # Test sorting in EventQueue
        # Mock SubVolume list
        subvolumes = [
            SubVolume((0, 0), 0.0, [0.0], 0.0),
            SubVolume((0, 1), 0.0, [0.0], 1.0),
            SubVolume((0, 2), 0.0, [0.0], 2.0),
            SubVolume((0, 3), 0.0, [0.0], 3.0),
            SubVolume((0, 4), 0.0, [0.0], 4.0),
        ]
        shuffle(subvolumes)

        # Sort SubVolume list
        event_queue: EventQueue = EventQueue(subvolumes)
        sorted_subvolumes = []
        while len(event_queue.heap) > 0:
            sorted_subvolumes.append(event_queue.extract_min())

        # Assert that the SubVolume list is sorted
        for i in range(len(sorted_subvolumes) - 1):
            self.assertTrue(sorted_subvolumes[i].next_event_time <= sorted_subvolumes[i + 1].next_event_time)
