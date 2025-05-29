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

def is_valid_move(current_state, next_state):
    # Only one box should move at a time
    differences = sum(1 for box in current_state if current_state[box] != next_state[box])
    if differences != 1:
        return False
    
    # The move should be to an adjacent cell
    for box in current_state:
        if current_state[box] != next_state[box]:
            if next_state[box] not in adjacency[current_state[box]]:
                return False
            
    # Check if target cell is occupied
    occupied_cells = set(next_state.values())
    if len(occupied_cells) != len(next_state):
        return False
            
    return True

def get_possible_moves(state):
    moves = []
    for box in state:
        # If box is at goal, don't move it
        if state[box] == goal_state[box]:
            continue
        
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            # Check if next_pos is occupied
            if next_pos in state.values():
                continue
                
            new_state = copy.deepcopy(state)
            new_state[box] = next_pos
            moves.append(new_state)
    
    return moves

def find_solution():
    queue = deque([(initial_state, [initial_state])])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, path = queue.popleft()
        
        if current_state == goal_state:
            return path
        
        for next_state in get_possible_moves(current_state):
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