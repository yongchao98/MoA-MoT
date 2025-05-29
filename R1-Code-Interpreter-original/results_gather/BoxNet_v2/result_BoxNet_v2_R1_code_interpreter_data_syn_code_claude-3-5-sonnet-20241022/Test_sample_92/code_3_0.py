from queue import PriorityQueue
import json
import copy

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

def get_next_states(current_state):
    next_states = []
    occupied_cells = set(current_state.values())
    
    # Move one box at a time
    for box in ["box1", "box2", "box3", "box4"]:
        # Skip if box is at goal
        if current_state[box] == goal_state[box]:
            continue
            
        current_pos = current_state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in occupied_cells:
                new_state = copy.deepcopy(current_state)
                new_state[box] = next_pos
                next_states.append(new_state)
    
    return next_states

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def solve():
    queue = [(initial_state, [initial_state])]
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, path = queue.pop(0)
        
        if current_state == goal_state:
            return path
            
        for next_state in get_next_states(current_state):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, path + [next_state]))
                
                # Early success check
                if next_state == goal_state:
                    return path + [next_state]
    
    return None

# Solve and print result
solution = solve()
if solution:
    # Ensure the solution is complete
    result = []
    for state in solution:
        if len(result) >= 30:  # Limit solution length to avoid overflow
            break
        result.append(state)
    print(json.dumps(result))
else:
    print("No solution found")