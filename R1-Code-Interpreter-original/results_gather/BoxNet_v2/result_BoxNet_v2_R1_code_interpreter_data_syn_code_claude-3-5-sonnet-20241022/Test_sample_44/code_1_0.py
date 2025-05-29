import heapq
import json
from copy import deepcopy

# Initial and goal states
initial_state = {
    "box1": "C1,2",
    "box2": "C2,4",
    "box3": "C2,3",
    "box4": "C1,4"
}

goal_state = {
    "box1": "C3,1",
    "box2": "C2,3",
    "box3": "C3,2",
    "box4": "C1,5"
}

# Adjacency dictionary
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"],
    "C2,6": ["C2,5", "C1,6", "C3,6"],
    "C3,1": ["C3,2", "C2,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5"],
    "C3,6": ["C3,5", "C2,6"]
}

class PriorityEntry:
    def __init__(self, priority, counter, item):
        self.priority = priority
        self.counter = counter
        self.item = item

    def __lt__(self, other):
        return (self.priority, self.counter) < (other.priority, other.counter)

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def heuristic(state):
    total = 0
    for box, pos in state.items():
        if pos != goal_state[box]:
            total += manhattan_distance(pos, goal_state[box])
    return total

def get_next_states(current_state):
    next_states = []
    occupied_cells = set(current_state.values())
    
    for box, pos in current_state.items():
        if pos == goal_state[box]:
            continue
            
        for next_pos in adjacency[pos]:
            if next_pos not in occupied_cells:
                new_state = deepcopy(current_state)
                new_state[box] = next_pos
                next_states.append(new_state)
    
    return next_states

def find_path():
    start_state = initial_state
    visited = set()
    counter = 0
    pq = [PriorityEntry(heuristic(start_state), counter, [start_state])]
    heapq.heapify(pq)
    
    while pq:
        entry = heapq.heappop(pq)
        path = entry.item
        current_state = path[-1]
        
        if current_state == goal_state:
            return path
            
        current_string = state_to_string(current_state)
        if current_string in visited:
            continue
            
        visited.add(current_string)
        
        for next_state in get_next_states(current_state):
            if state_to_string(next_state) not in visited:
                counter += 1
                new_path = path + [next_state]
                new_cost = len(new_path) - 1
                new_priority = new_cost + heuristic(next_state)
                heapq.heappush(pq, PriorityEntry(new_priority, counter, new_path))
    
    return None

# Find and print solution
solution = find_path()
if solution:
    print("<<<" + json.dumps(solution) + ">>>")
else:
    print("No solution found")