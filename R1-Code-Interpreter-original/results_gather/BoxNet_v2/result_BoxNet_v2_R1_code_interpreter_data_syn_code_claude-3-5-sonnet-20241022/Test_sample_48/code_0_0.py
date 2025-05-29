import json
from collections import deque
import copy

# Initial and goal states
initial_state = {
    "box1": "C5,5",
    "box2": "C1,3",
    "box3": "C3,3",
    "box4": "C4,5",
    "box5": "C4,4"
}

goal_state = {
    "box1": "C2,5",
    "box2": "C3,4",
    "box3": "C1,5",
    "box4": "C3,1",
    "box5": "C1,3"
}

# Adjacency dictionary
adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C2,5"], "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"], "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"], "C2,5": ["C2,4", "C1,5", "C3,5"],
    "C3,1": ["C3,2", "C2,1", "C4,1"], "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C2,5", "C4,5"], "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"], "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"], "C4,5": ["C4,4", "C3,5", "C5,5"],
    "C5,1": ["C5,2", "C4,1"], "C5,2": ["C5,1", "C5,3", "C4,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3"], "C5,4": ["C5,3", "C5,5", "C4,4"],
    "C5,5": ["C5,4", "C4,5"]
}

def state_to_string(state):
    return str(sorted(state.items()))

def is_valid_move(current_state, next_state):
    # Check if only one box moved to an adjacent cell
    differences = 0
    moved_box = None
    for box in current_state:
        if current_state[box] != next_state[box]:
            differences += 1
            moved_box = box
            if differences > 1:
                return False
            if next_state[box] not in adjacency[current_state[box]]:
                return False
    
    # Check if destination cell is occupied by another box
    for box in next_state:
        if box != moved_box and next_state[box] == next_state[moved_box]:
            return False
    
    return True

def get_next_states(current_state):
    next_states = []
    for box in current_state:
        # Skip if box is already at goal
        if current_state[box] == goal_state[box]:
            continue
        
        # Try moving the box to adjacent cells
        for adj_cell in adjacency[current_state[box]]:
            new_state = copy.deepcopy(current_state)
            new_state[box] = adj_cell
            
            # Check if move is valid
            if is_valid_move(current_state, new_state):
                # Check if destination is not occupied
                if adj_cell not in new_state.values() or new_state[box] == adj_cell:
                    next_states.append(new_state)
    
    return next_states

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def heuristic(state):
    total_distance = 0
    for box in state:
        total_distance += manhattan_distance(state[box], goal_state[box])
    return total_distance

def find_solution():
    start_state = initial_state
    visited = set()
    queue = deque([(start_state, [start_state])])
    visited.add(state_to_string(start_state))
    
    while queue:
        current_state, path = queue.popleft()
        
        if current_state == goal_state:
            return path
        
        next_states = get_next_states(current_state)
        # Sort next states by heuristic value
        next_states.sort(key=heuristic)
        
        for next_state in next_states:
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, path + [next_state]))
    
    return None

# Find and print solution
solution = find_solution()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")