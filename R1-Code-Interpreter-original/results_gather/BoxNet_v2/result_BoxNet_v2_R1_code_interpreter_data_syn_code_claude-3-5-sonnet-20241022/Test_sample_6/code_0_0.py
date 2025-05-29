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
    for box in current_state:
        if current_state[box] != next_state[box]:
            if next_state[box] not in adjacency[current_state[box]]:
                return False
            differences += 1
    return differences <= 1

def get_manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def get_total_distance(state):
    total = 0
    for box in state:
        if state[box] != goal_state[box]:
            total += get_manhattan_distance(state[box], goal_state[box])
    return total

def get_next_states(current_state):
    next_states = []
    for box in current_state:
        if current_state[box] == goal_state[box]:
            continue
        current_pos = current_state[box]
        for next_pos in adjacency[current_pos]:
            # Check if next_pos is occupied by another box
            if next_pos not in current_state.values():
                new_state = copy.deepcopy(current_state)
                new_state[box] = next_pos
                next_states.append(new_state)
    return next_states

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
        next_states.sort(key=get_total_distance)  # Heuristic ordering
        
        for next_state in next_states:
            state_str = state_to_string(next_state)
            if state_str not in visited and is_valid_move(current_state, next_state):
                visited.add(state_str)
                queue.append((next_state, path + [next_state]))
    
    return None

# Find and print solution
solution = find_solution()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")