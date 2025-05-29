import heapq
from copy import deepcopy
import json

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
    def __init__(self, priority, counter, data):
        self.priority = priority
        self.counter = counter
        self.data = data

    def __lt__(self, other):
        return (self.priority, self.counter) < (other.priority, other.counter)

def manhattan_distance(pos1, pos2):
    row1, col1 = map(int, pos1.replace('C', '').split(','))
    row2, col2 = map(int, pos2.replace('C', '').split(','))
    return abs(row1 - row2) + abs(col1 - col2)

def heuristic(state):
    total = 0
    for box, pos in state.items():
        if pos != goal_state[box]:
            total += manhattan_distance(pos, goal_state[box])
    return total

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

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
    counter = 0
    start_h = heuristic(initial_state)
    pq = [PriorityEntry(start_h, counter, [initial_state])]
    visited = {state_to_string(initial_state)}
    
    while pq:
        entry = heapq.heappop(pq)
        cost = entry.priority - heuristic(entry.data[-1])
        path = entry.data
        current_state = path[-1]
        
        if all(current_state[box] == goal_state[box] for box in current_state):
            return path
            
        for next_state in get_next_states(current_state):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                h = heuristic(next_state)
                new_path = path + [next_state]
                counter += 1
                heapq.heappush(pq, PriorityEntry(cost + 1 + h, counter, new_path))
    
    return None

# Find and print the solution
solution = find_path()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")