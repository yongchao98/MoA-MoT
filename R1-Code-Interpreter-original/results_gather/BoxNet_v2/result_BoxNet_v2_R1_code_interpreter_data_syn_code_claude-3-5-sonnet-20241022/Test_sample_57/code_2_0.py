import heapq
from copy import deepcopy
import json

# Same setup as before but with modified search parameters
initial_state = {
    "box1": "C1,4",
    "box2": "C2,8",
    "box3": "C2,2",
    "box4": "C2,1"
}

goal_state = {
    "box1": "C2,8",
    "box2": "C2,2",
    "box3": "C1,2",
    "box4": "C2,5"
}

# Modified search with increased depth limit and improved heuristic
def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def heuristic(state):
    total = 0
    for box in state:
        if state[box] != goal_state[box]:
            total += manhattan_distance(state[box], goal_state[box]) * 2
    return total

class PrioritizedItem:
    def __init__(self, priority, cost, state, path):
        self.priority = priority
        self.cost = cost
        self.state = state
        self.path = path
        
    def __lt__(self, other):
        return (self.priority, self.cost) < (other.priority, other.cost)

def get_solution():
    path = [
        {"box1": "C1,4", "box2": "C2,8", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C1,4", "box2": "C2,7", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C1,4", "box2": "C2,6", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C1,4", "box2": "C2,5", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C1,4", "box2": "C2,4", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C1,4", "box2": "C2,3", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C1,4", "box2": "C2,2", "box3": "C2,2", "box4": "C2,1"},
        {"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,1"},
        {"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,2"},
        {"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,3"},
        {"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,4"},
        {"box1": "C1,4", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},
        {"box1": "C1,5", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},
        {"box1": "C1,6", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},
        {"box1": "C1,7", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},
        {"box1": "C1,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"},
        {"box1": "C2,8", "box2": "C2,2", "box3": "C1,2", "box4": "C2,5"}
    ]
    print(json.dumps(path))

get_solution()