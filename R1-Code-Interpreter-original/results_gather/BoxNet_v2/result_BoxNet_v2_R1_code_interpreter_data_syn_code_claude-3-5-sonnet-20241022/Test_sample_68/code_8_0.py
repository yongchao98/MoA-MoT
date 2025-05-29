import json
from collections import deque

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
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"], "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1"], "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"], "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"], "C2,6": ["C2,5", "C2,7", "C1,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7"], "C2,8": ["C2,7", "C1,8"]
}

def get_next_states(state):
    next_states = []
    # Try moving each box that's not at its goal
    for box in state:
        if state[box] == goal_state[box]:
            continue
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if next_pos not in state.values():
                new_state = state.copy()
                new_state[box] = next_pos
                next_states.append(new_state)
    return next_states

def solve():
    start = initial_state
    queue = deque([(start, [start])])
    visited = {json.dumps(start, sort_keys=True)}
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if current state is goal state
        if current_state == goal_state:
            return path
        
        # Try all possible next states
        for next_state in get_next_states(current_state):
            state_str = json.dumps(next_state, sort_keys=True)
            if state_str not in visited:
                visited.add(state_str)
                new_path = path + [next_state]
                queue.append((next_state, new_path))
                
                # Early success check
                if next_state == goal_state:
                    return new_path
    
    return None

# Get and format solution
solution = solve()
if solution:
    # Format each state in the solution
    formatted_solution = []
    for state in solution:
        # Ensure consistent ordering of items in each state
        formatted_state = dict(sorted(state.items()))
        formatted_solution.append(formatted_state)
    
    # Print the complete solution
    print(json.dumps(formatted_solution))
else:
    print("No solution found")