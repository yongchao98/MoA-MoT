import json
from collections import deque
import copy

# Problem setup
initial_state = {
    "box1": "C1,3",
    "box2": "C2,4",
    "box3": "C1,6",
    "box4": "C2,2",
    "box5": "C1,4"
}

goal_state = {
    "box1": "C2,1",
    "box2": "C1,2",
    "box3": "C1,5",
    "box4": "C2,2",
    "box5": "C2,6"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"],
    "C2,6": ["C2,5", "C1,6"]
}

def state_to_string(state):
    return str(sorted(state.items()))

def is_goal_reached(state):
    return all(state[box] == goal_state[box] for box in state)

def get_valid_moves(state):
    moves = []
    occupied_cells = set(state.values())
    
    for box, current_pos in state.items():
        # If box is at goal, skip it
        if current_pos == goal_state[box]:
            continue
            
        # Check each adjacent cell
        for next_pos in adjacency[current_pos]:
            # Check if the adjacent cell is not occupied
            if next_pos not in occupied_cells:
                new_state = copy.deepcopy(state)
                new_state[box] = next_pos
                moves.append(new_state)
    
    return moves

def manhattan_distance(pos1, pos2):
    row1, col1 = int(pos1.split(',')[0][1]), int(pos1.split(',')[1])
    row2, col2 = int(pos2.split(',')[0][1]), int(pos2.split(',')[1])
    return abs(row1 - row2) + abs(col1 - col2)

def heuristic(state):
    return sum(manhattan_distance(state[box], goal_state[box]) for box in state)

def find_solution():
    start_state = initial_state
    visited = {state_to_string(start_state)}
    queue = deque([(start_state, [start_state])])
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal_reached(current_state):
            return path
            
        for next_state in get_valid_moves(current_state):
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