from collections import deque

# Initial and goal states
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

# Adjacency list
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

# Function to generate next states
def generate_next_states(state):
    next_states = []
    for box, location in state.items():
        if location == goal_state[box]:
            continue
        for neighbor in adjacency[location]:
            new_state = state.copy()
            new_state[box] = neighbor
            if len(set(new_state.values())) == len(new_state):  # Ensure no two boxes in the same cell
                next_states.append(new_state)
    return next_states

# BFS to find the path
def bfs(initial_state, goal_state):
    queue = deque([(initial_state, [initial_state])])
    visited = set()
    visited.add(tuple(initial_state.items()))
    
    while queue:
        current_state, path = queue.popleft()
        
        if current_state == goal_state:
            return path
        
        for next_state in generate_next_states(current_state):
            state_tuple = tuple(next_state.items())
            if state_tuple not in visited:
                visited.add(state_tuple)
                queue.append((next_state, path + [next_state]))
    
    return None

# Find the plan
plan = bfs(initial_state, goal_state)
print(plan)