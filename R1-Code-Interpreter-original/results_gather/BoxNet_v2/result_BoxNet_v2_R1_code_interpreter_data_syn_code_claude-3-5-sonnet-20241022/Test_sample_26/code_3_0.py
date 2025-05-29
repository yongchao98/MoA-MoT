import json
from copy import deepcopy

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
    occupied = set(current_state.values())
    
    # Try moving each box that hasn't reached its goal
    for box in current_state:
        if current_state[box] == goal_state[box]:
            continue
        
        curr_pos = current_state[box]
        for next_pos in adjacency[curr_pos]:
            if next_pos not in occupied:
                new_state = deepcopy(current_state)
                new_state[box] = next_pos
                next_states.append(new_state)
    
    return next_states

def solve():
    # Start with initial state
    states = [initial_state]
    visited = {json.dumps(initial_state, sort_keys=True)}
    queue = [(states, 0)]
    
    while queue:
        current_states, depth = queue.pop(0)
        current = current_states[-1]
        
        # Check if we've reached the goal
        if current == goal_state:
            return current_states
            
        # Limit search depth
        if depth > 20:
            continue
            
        # Get next possible states
        for next_state in get_next_states(current):
            state_str = json.dumps(next_state, sort_keys=True)
            if state_str not in visited:
                visited.add(state_str)
                new_states = current_states + [next_state]
                queue.append((new_states, depth + 1))
        
        # Keep queue size manageable
        if len(queue) > 500:
            queue = queue[:500]
    
    return None

# Solve and print result
solution = solve()
if solution:
    print(f"<<<{json.dumps(solution)}>>>")
else:
    print("No solution found")