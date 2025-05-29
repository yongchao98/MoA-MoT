from enum import Enum
from typing import Tuple, Dict

class Orientation(Enum):
    north = 0
    west = 1
    south = 2
    east = 3
    half = 4
    empty = 5
    blocked = 6

Tatami = Tuple[Orientation, int]

def orientation(tatami: Tatami) -> Orientation:
    return tatami[0]

class Room:
    def __init__(self, height: int, width: int) -> None:
        self.width = width
        self.height = height
        self.tiles = [[(Orientation.empty, -1) for j in range(width)] for i in range(height)]
        self.corners = [[0 for j in range(width+1)] for i in range(height+1)]

    def is_empty_spot(self, pos: Tuple[int, int]) -> bool:
        return self.tiles[pos[0]][pos[1]][0] == Orientation.empty

    def can_place_tatami(self, pos: Tuple[int, int], tatami: Tatami) -> bool:
        if not self.is_empty_spot(pos):
            return False
        corners: Dict[str, int] = self.number_of_corners(pos)
        if orientation(tatami) == Orientation.half and (corners["se"] > 2 or corners["ne"] > 2 or corners["nw"] > 2 or corners["sw"] > 2):
            return False
        return True

    def number_of_corners(self, pos: Tuple[int, int]) -> Dict[str, int]:
        corners: Dict[str, int] = {
            "nw": self.corners[pos[0]+1][pos[1]],
            "ne": self.corners[pos[0]+1][pos[1]+1],
            "sw": self.corners[pos[0]][pos[1]],
            "se": self.corners[pos[0]][pos[1]+1]
        }
        return corners

def main_solution(height: int, width: int, tatami_type: str, tatami_index: int) -> bool:
    room = Room(height, width)
    tatami = (Orientation[tatami_type], tatami_index)
    for i in range(height):
        for j in range(width):
            if room.can_place_tatami((i, j), tatami):
                return True
    return False

# Test the function with the given input
result = main_solution(6, 10, 'half', 2)
print(result)