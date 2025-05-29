from collections import deque

# Initial and goal states
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

# Adjacency list
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

# Function to generate the plan
def generate_plan(initial_state, goal_state, adjacency):
    queue = deque([(initial_state, [])])
    visited = set()
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if the current state is the goal state
        if current_state == goal_state:
            return path + [current_state]
        
        # Generate possible moves
        for box, current_location in current_state.items():
            if current_location != goal_state[box]:  # Only move if not at goal
                for neighbor in adjacency[current_location]:
                    new_state = current_state.copy()
                    new_state[box] = neighbor
                    state_tuple = tuple(sorted(new_state.items()))
                    
                    if state_tuple not in visited:
                        visited.add(state_tuple)
                        queue.append((new_state, path + [current_state]))
    
    return []

# Generate the plan
plan = generate_plan(initial_state, goal_state, adjacency)
print(plan)