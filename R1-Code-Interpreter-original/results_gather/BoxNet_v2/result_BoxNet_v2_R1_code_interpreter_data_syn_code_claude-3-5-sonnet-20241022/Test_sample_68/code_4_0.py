import json
from collections import deque

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

def is_valid_move(state, box, new_pos):
    # Check if new position is already occupied
    return new_pos not in state.values()

def get_next_states(state):
    next_states = []
    # Try moving each box that's not at its goal
    for box in state:
        if state[box] == goal_state[box]:
            continue
        current_pos = state[box]
        for next_pos in adjacency[current_pos]:
            if is_valid_move(state, box, next_pos):
                new_state = state.copy()
                new_state[box] = next_pos
                next_states.append(new_state)
    return next_states

def state_to_string(state):
    return ','.join(f"{k}:{v}" for k, v in sorted(state.items()))

def solve():
    queue = deque([(initial_state, [initial_state])])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if current state is goal state
        if all(current_state[box] == goal_state[box] for box in current_state):
            return path
        
        # Get and process all possible next states
        for next_state in get_next_states(current_state):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                new_path = path + [next_state]
                queue.append((next_state, new_path))
    
    return None

# Find solution and format output
solution = solve()
if solution:
    # Ensure the solution is properly formatted
    formatted_solution = []
    for state in solution:
        formatted_state = dict(sorted(state.items()))
        formatted_solution.append(formatted_state)
    
    # Print the complete solution
    result = json.dumps(formatted_solution, separators=(',', ':'))
    print(result)
else:
    print("No solution found")