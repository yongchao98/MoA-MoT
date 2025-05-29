from collections import deque

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

# Function to generate the plan
def generate_plan(initial_state, goal_state, adjacency):
    queue = deque([(initial_state, [])])
    visited = set()
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if the current state is the goal state
        if current_state == goal_state:
            return path + [current_state]
        
        # Generate all possible moves
        for box, current_position in current_state.items():
            if current_position == goal_state[box]:
                continue  # Skip if the box is already at its goal
            
            for neighbor in adjacency[current_position]:
                # Create a new state with the box moved to the neighbor
                new_state = current_state.copy()
                new_state[box] = neighbor
                
                # Check if the new state is valid (no two boxes in the same cell)
                if len(set(new_state.values())) == len(new_state):
                    state_tuple = tuple(sorted(new_state.items()))
                    if state_tuple not in visited:
                        visited.add(state_tuple)
                        queue.append((new_state, path + [current_state]))
    
    return []

# Generate the plan
plan = generate_plan(initial_state, goal_state, adjacency)
print(plan)