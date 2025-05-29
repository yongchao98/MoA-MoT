from collections import deque
import json
import copy

# Initial and goal states
initial_state = {
    "box1": "C1,2",
    "box2": "C2,4",
    "box3": "C2,3",
    "box4": "C1,4"
}

goal_state = {
    "box1": "C3,1",
    "box2": "C2,3",
    "box3": "C3,2",
    "box4": "C1,5"
}

# Adjacency dictionary
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"],
    "C2,6": ["C2,5", "C1,6", "C3,6"],
    "C3,1": ["C3,2", "C2,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5"],
    "C3,6": ["C3,5", "C2,6"]
}

def is_valid_state(state):
    # Check if any two boxes share the same cell
    positions = state.values()
    return len(positions) == len(set(positions))

def get_next_states(current_state):
    next_states = []
    
    # Try moving each box
    for box in current_state:
        # If box is at goal, skip it
        if current_state[box] == goal_state[box]:
            continue
            
        current_pos = current_state[box]
        # Try each adjacent position
        for next_pos in adjacency[current_pos]:
            new_state = copy.deepcopy(current_state)
            new_state[box] = next_pos
            if is_valid_state(new_state):
                next_states.append(new_state)
    
    return next_states

def state_to_string(state):
    return json.dumps(state, sort_keys=True)

def find_solution():
    start_state = initial_state
    visited = set()
    queue = deque([(start_state, [start_state])])
    visited.add(state_to_string(start_state))
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if we reached the goal
        all_at_goal = all(current_state[box] == goal_state[box] for box in current_state)
        if all_at_goal:
            return path
            
        # Generate and check next states
        for next_state in get_next_states(current_state):
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