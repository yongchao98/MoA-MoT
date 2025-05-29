from collections import deque

# Initial and goal states
initial_state = {
    "box1": "C4,3",
    "box2": "C4,6",
    "box3": "C2,4",
    "box4": "C1,2"
}

goal_state = {
    "box1": "C2,1",
    "box2": "C1,4",
    "box3": "C1,6",
    "box4": "C4,2"
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
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"],
    "C3,6": ["C3,5", "C2,6", "C4,6"],
    "C4,1": ["C4,2", "C3,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5"],
    "C4,6": ["C4,5", "C3,6"]
}

# Function to generate the plan
def generate_plan(initial_state, goal_state, adjacency):
    queue = deque([(initial_state, [])])  # (current_state, path)
    visited = set()
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if the current state is the goal state
        if current_state == goal_state:
            return path + [current_state]
        
        # Mark the current state as visited
        state_tuple = tuple(sorted(current_state.items()))
        if state_tuple in visited:
            continue
        visited.add(state_tuple)
        
        # Generate new states by moving each box
        for box, current_cell in current_state.items():
            if current_cell == goal_state[box]:
                continue  # Skip if the box is already at its goal
            
            for neighbor in adjacency[current_cell]:
                # Check if the neighbor cell is free
                if neighbor not in current_state.values():
                    new_state = current_state.copy()
                    new_state[box] = neighbor
                    queue.append((new_state, path + [current_state]))
    
    return []  # Return an empty list if no plan is found

# Generate the plan
plan = generate_plan(initial_state, goal_state, adjacency)
print(plan)