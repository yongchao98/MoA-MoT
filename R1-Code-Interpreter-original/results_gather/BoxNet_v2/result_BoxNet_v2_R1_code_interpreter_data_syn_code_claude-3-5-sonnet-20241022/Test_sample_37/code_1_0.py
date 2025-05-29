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

class PrioritizedState:
    def __init__(self, priority, steps, state, path):
        self.priority = priority
        self.steps = steps
        self.state = state
        self.path = path
        
    def __lt__(self, other):
        return (self.priority, self.steps) < (other.priority, other.steps)

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

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

def get_valid_moves(state):
    moves = []
    occupied_cells = set(state.values())
    
    for box, current_pos in state.items():
        if current_pos == goal_state[box]:
            continue
            
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied_cells:
                new_state = deepcopy(state)
                new_state[box] = next_pos
                moves.append(new_state)
    
    return moves

def find_path():
    start_string = state_to_string(initial_state)
    visited = {start_string}
    pq = [PrioritizedState(heuristic(initial_state), 0, initial_state, [initial_state])]
    heapq.heapify(pq)
    
    while pq:
        current = heapq.heappop(pq)
        current_state = current.state
        
        if current_state == goal_state:
            return current.path
            
        for next_state in get_valid_moves(current_state):
            state_string = state_to_string(next_state)
            if state_string not in visited:
                visited.add(state_string)
                new_path = current.path + [next_state]
                priority = len(new_path) + heuristic(next_state)
                heapq.heappush(pq, PrioritizedState(priority, len(new_path), next_state, new_path))
    
    return None

# Find and print solution
solution = find_path()
if solution:
    print("<<<" + json.dumps(solution) + ">>>")
else:
    print("No solution found")