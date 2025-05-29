import json
from copy import deepcopy
from heapq import heappush, heappop

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1[1:].split(','))
    x2, y2 = map(int, pos2[1:].split(','))
    return abs(x1 - x2) + abs(y1 - y2)

class PrioritizedState:
    def __init__(self, priority, steps, state, path):
        self.priority = priority
        self.steps = steps
        self.state = state
        self.path = path
        
    def __lt__(self, other):
        return (self.priority, self.steps) < (other.priority, other.steps)

initial_state = {
    "box1": "C5,5", "box2": "C1,3", "box3": "C3,3",
    "box4": "C4,5", "box5": "C4,4"
}

goal_state = {
    "box1": "C2,5", "box2": "C3,4", "box3": "C1,5",
    "box4": "C3,1", "box5": "C1,3"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C2,5"], "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C1,5", "C3,5"], "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C2,5", "C4,5"], "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C3,5", "C5,5"], "C5,1": ["C5,2", "C4,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2"], "C5,3": ["C5,2", "C5,4", "C4,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4"], "C5,5": ["C5,4", "C4,5"]
}

def get_valid_moves(state):
    moves = []
    occupied = set(state.values())
    
    for box, pos in state.items():
        if pos == goal_state[box]:
            continue
            
        for next_pos in adjacency[pos]:
            if next_pos not in occupied:
                new_state = deepcopy(state)
                new_state[box] = next_pos
                moves.append(new_state)
    
    return moves

def is_goal_reached(state):
    return all(state[box] == goal_state[box] for box in state)

def find_solution():
    max_steps = 50  # Limit the number of steps to prevent infinite loops
    visited = set()
    queue = [(0, 0, initial_state, [initial_state])]  # (steps, id, state, path)
    
    while queue:
        steps, _, current_state, path = heappop(queue)
        
        if steps > max_steps:
            continue
            
        if is_goal_reached(current_state):
            return path
            
        state_hash = tuple(sorted(current_state.items()))
        if state_hash in visited:
            continue
            
        visited.add(state_hash)
        
        for next_state in get_valid_moves(current_state):
            new_path = path + [next_state]
            heappush(queue, (steps + 1, len(queue), next_state, new_path))
    
    return None

# Find and print solution
solution = find_solution()
if solution:
    # Format the solution as a JSON string with proper indentation
    formatted_solution = json.dumps(solution, indent=None)
    print(formatted_solution)
else:
    print("No solution found")