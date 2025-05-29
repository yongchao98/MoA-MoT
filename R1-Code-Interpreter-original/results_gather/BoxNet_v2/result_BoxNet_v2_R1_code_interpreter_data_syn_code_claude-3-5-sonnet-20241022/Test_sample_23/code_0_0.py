import copy
from collections import deque

def get_valid_moves(current_state, goal_state, adjacency):
    valid_moves = []
    # For each box
    for box in current_state:
        # If box is already at goal, skip it
        if current_state[box] == goal_state[box]:
            continue
        
        current_pos = current_state[box]
        # Check each adjacent cell
        for next_pos in adjacency[current_pos]:
            # Check if the adjacent cell is occupied
            if next_pos not in current_state.values():
                # Create new state with this move
                new_state = copy.deepcopy(current_state)
                new_state[box] = next_pos
                valid_moves.append(new_state)
    
    return valid_moves

def state_to_string(state):
    return str(sorted(state.items()))

def find_solution(initial_state, goal_state, adjacency):
    queue = deque([(initial_state, [initial_state])])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if we reached the goal
        if all(current_state[box] == goal_state[box] for box in current_state):
            return path
        
        # Get valid moves from current state
        for next_state in get_valid_moves(current_state, goal_state, adjacency):
            state_str = state_to_string(next_state)
            if state_str not in visited:
                visited.add(state_str)
                queue.append((next_state, path + [next_state]))
    
    return None

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

# Adjacency list
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

# Find solution
solution = find_solution(initial_state, goal_state, adjacency)

# Print solution in required format
if solution:
    print(str(solution).replace("'", '"'))