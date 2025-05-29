import json
from collections import deque
from copy import deepcopy

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

def get_next_states(state):
    next_states = []
    occupied = set(state.values())
    
    for box in state:
        current_pos = state[box]
        
        # If box is at goal, skip it
        if current_pos == goal_state[box]:
            continue
            
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied or next_pos == current_pos:
                new_state = state.copy()
                new_state[box] = next_pos
                next_states.append(new_state)
    
    return next_states

def state_to_string(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))

def is_valid_transition(state1, state2):
    diff_count = sum(1 for k in state1 if state1[k] != state2[k])
    if diff_count != 1:
        return False
    
    # Check if the moved box went to an adjacent cell
    for k in state1:
        if state1[k] != state2[k]:
            return state2[k] in adjacency[state1[k]]
    
    return True

def bfs_with_path():
    queue = deque([(initial_state, [initial_state])])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, path = queue.popleft()
        
        if current_state == goal_state:
            return path
        
        for next_state in get_next_states(current_state):
            state_str = state_to_string(next_state)
            
            if state_str not in visited and is_valid_transition(current_state, next_state):
                visited.add(state_str)
                new_path = path + [next_state]
                queue.append((next_state, new_path))
                
                # Early success check
                if next_state == goal_state:
                    return new_path
        
        # Limit search space
        if len(path) > 50:
            continue
    
    return None

def verify_solution(solution):
    if not solution:
        return False
    
    # Check first and last states
    if solution[0] != initial_state or solution[-1] != goal_state:
        return False
    
    # Check transitions
    for i in range(len(solution) - 1):
        if not is_valid_transition(solution[i], solution[i + 1]):
            return False
    
    return True

# Find solution
solution = bfs_with_path()

# Format and verify solution
if solution and verify_solution(solution):
    formatted_solution = []
    for state in solution:
        formatted_solution.append(dict(sorted(state.items())))
    print(json.dumps(formatted_solution))
else:
    print("No valid solution found")