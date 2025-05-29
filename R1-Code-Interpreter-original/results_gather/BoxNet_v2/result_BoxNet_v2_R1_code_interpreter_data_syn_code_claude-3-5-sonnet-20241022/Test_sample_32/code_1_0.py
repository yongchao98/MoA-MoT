import heapq
import json
from copy import deepcopy

# Initial and goal states
initial_state = {
    "box1": "C4,3",
    "box2": "C4,6",
    "box3": "C2,4",
    "box4": "C1,2"
}

goal_state = {
    "box1": "C2,1",
    "box2": "C1,4",
    "box3": "C1,6",
    "box4": "C4,2"
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
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"],
    "C3,6": ["C3,5", "C2,6", "C4,6"],
    "C4,1": ["C4,2", "C3,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5"],
    "C4,6": ["C4,5", "C3,6"]
}

class PriorityEntry:
    def __init__(self, priority, counter, path, state):
        self.priority = priority
        self.counter = counter
        self.path = path
        self.state = state

    def __lt__(self, other):
        return (self.priority, self.counter) < (other.priority, other.counter)

def manhattan_distance(pos1, pos2):
    row1, col1 = map(int, pos1[1:].split(','))
    row2, col2 = map(int, pos2[1:].split(','))
    return abs(row1 - row2) + abs(col1 - col2)

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def get_valid_moves(current_state):
    valid_moves = []
    occupied_cells = set(current_state.values())
    
    for box, pos in current_state.items():
        # If box is at goal, skip it
        if pos == goal_state[box]:
            continue
            
        # Check each adjacent cell
        for next_pos in adjacency[pos]:
            if next_pos not in occupied_cells:
                new_state = deepcopy(current_state)
                new_state[box] = next_pos
                valid_moves.append(new_state)
    
    return valid_moves

def heuristic(state):
    total = 0
    for box, pos in state.items():
        if pos != goal_state[box]:
            total += manhattan_distance(pos, goal_state[box])
    return total

def find_solution():
    start_state = initial_state
    visited = set()
    counter = 0
    pq = [PriorityEntry(heuristic(start_state), counter, [start_state], start_state)]
    heapq.heapify(pq)
    
    while pq:
        current = heapq.heappop(pq)
        current_state = current.state
        path = current.path
        
        state_str = state_to_string(current_state)
        
        if current_state == goal_state:
            return path
            
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        for next_state in get_valid_moves(current_state):
            if state_to_string(next_state) not in visited:
                counter += 1
                new_path = path + [next_state]
                new_cost = len(new_path) - 1
                priority = new_cost + heuristic(next_state)
                heapq.heappush(pq, PriorityEntry(priority, counter, new_path, next_state))
    
    return None

# Find and print solution
solution = find_solution()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")