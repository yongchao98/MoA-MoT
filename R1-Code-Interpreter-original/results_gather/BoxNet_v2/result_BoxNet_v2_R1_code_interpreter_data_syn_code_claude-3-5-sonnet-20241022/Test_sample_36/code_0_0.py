import json
from collections import deque
import copy

# Initial and goal states
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

# Adjacency dictionary
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

def state_to_string(state):
    return str(sorted(state.items()))

def is_valid_move(current_state, box, new_pos):
    # Check if the new position is occupied by another box
    for other_box, pos in current_state.items():
        if other_box != box and pos == new_pos:
            return False
    return True

def get_next_states(current_state):
    next_states = []
    for box in current_state:
        # If box is already at goal, skip it
        if current_state[box] == goal_state[box]:
            continue
        
        current_pos = current_state[box]
        for next_pos in adjacency[current_pos]:
            if is_valid_move(current_state, box, next_pos):
                new_state = copy.deepcopy(current_state)
                new_state[box] = next_pos
                next_states.append(new_state)
    return next_states

def manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

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
        
        if all(current_state[box] == goal_state[box] for box in current_state):
            return path
            
        next_states = get_next_states(current_state)
        # Sort next states by heuristic value
        next_states.sort(key=lambda x: heuristic(x))
        
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