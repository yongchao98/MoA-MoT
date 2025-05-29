import heapq
from copy import deepcopy
import json

# Problem setup
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
    # First try to move boxes that aren't at their goals
    for box in current_state:
        if current_state[box] == goal_state[box]:
            continue
        current_pos = current_state[box]
        for next_pos in adjacency[current_pos]:
            # Check if the next position is empty
            if next_pos not in current_state.values():
                new_state = deepcopy(current_state)
                new_state[box] = next_pos
                next_states.append(new_state)
    return next_states

def solve():
    from collections import deque
    
    def state_to_string(state):
        return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))
    
    queue = deque([(initial_state, [initial_state])])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, path = queue.popleft()
        
        if current_state == goal_state:
            return path
            
        # Generate and process next states
        for next_state in get_next_states(current_state):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                new_path = path + [next_state]
                queue.append((next_state, new_path))
                
                # Early success check
                if next_state == goal_state:
                    return new_path
    
    return None

# Solve and print result
solution = solve()
if solution:
    # Format solution to ensure it's complete
    result = []
    for state in solution:
        result.append(dict(sorted(state.items())))
    print(json.dumps(result))
else:
    print("No solution found")