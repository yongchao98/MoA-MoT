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

def get_next_states(current_state):
    next_states = []
    
    # For each box
    for box in current_state:
        # Skip if box is already at goal
        if current_state[box] == goal_state[box]:
            continue
            
        current_pos = current_state[box]
        # Try each adjacent position
        for next_pos in adjacency[current_pos]:
            # Check if next position is occupied
            if next_pos not in current_state.values():
                new_state = current_state.copy()
                new_state[box] = next_pos
                next_states.append(new_state)
                
    return next_states

def state_to_string(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))

def find_solution():
    queue = [(initial_state, [initial_state])]
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, path = queue.pop(0)
        
        # Check if we reached goal state
        if current_state == goal_state:
            return path
            
        # Generate next possible states
        next_states = get_next_states(current_state)
        
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