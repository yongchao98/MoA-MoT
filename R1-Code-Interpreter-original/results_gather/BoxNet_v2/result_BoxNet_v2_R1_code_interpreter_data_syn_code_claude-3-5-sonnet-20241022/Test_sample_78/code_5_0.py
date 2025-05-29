import heapq
from copy import deepcopy
import json

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

adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"],
    "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"],
    "C2,6": ["C2,5", "C2,7", "C1,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7"],
    "C2,8": ["C2,7", "C1,8"]
}

class PQEntry:
    def __init__(self, priority, moves, state_str, state, path):
        self.priority = priority
        self.moves = moves
        self.state_str = state_str
        self.state = state
        self.path = path

    def __lt__(self, other):
        return (self.priority, self.moves) < (other.priority, other.moves)

def manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def get_next_states(current_state):
    next_states = []
    occupied_cells = set(current_state.values())
    
    for box in current_state:
        if current_state[box] == goal_state[box]:
            continue
            
        current_pos = current_state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied_cells:
                new_state = deepcopy(current_state)
                new_state[box] = next_pos
                next_states.append(new_state)
    
    return next_states

def heuristic(state):
    total = 0
    for box in state:
        if state[box] != goal_state[box]:
            total += manhattan_distance(state[box], goal_state[box])
    return total

def solve():
    visited = set()
    start_state = initial_state
    start_str = state_to_string(start_state)
    h = heuristic(start_state)
    
    # Create initial entry
    start_entry = PQEntry(h, 0, start_str, start_state, [start_state])
    queue = [start_entry]
    visited.add(start_str)
    
    while queue:
        current = heapq.heappop(queue)
        current_state = current.state
        
        if current_state == goal_state:
            return current.path
        
        for next_state in get_next_states(current_state):
            next_str = state_to_string(next_state)
            if next_str not in visited:
                visited.add(next_str)
                new_path = current.path + [next_state]
                moves = len(new_path) - 1
                h = heuristic(next_state)
                priority = moves + h
                
                entry = PQEntry(priority, moves, next_str, next_state, new_path)
                heapq.heappush(queue, entry)
    
    return None

# Solve and print result
solution = solve()
if solution:
    print("<<<" + json.dumps(solution) + ">>>")
else:
    print("No solution found")