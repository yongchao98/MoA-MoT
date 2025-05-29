from collections import deque

# Define the initial and goal states
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

# Define the adjacency list
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

# Function to perform BFS to find the shortest path for a box
def bfs_path(start, goal, current_state):
    queue = deque([(start, [start])])
    visited = set()
    
    while queue:
        position, path = queue.popleft()
        if position == goal:
            return path
        if position not in visited:
            visited.add(position)
            for neighbor in adjacency[position]:
                if neighbor not in current_state.values():  # Ensure no box is in the neighbor
                    queue.append((neighbor, path + [neighbor]))
    return []

# Plan generation
def generate_plan(initial_state, goal_state):
    current_state = initial_state.copy()
    plan = [current_state.copy()]
    
    for box, goal_position in goal_state.items():
        if current_state[box] != goal_position:
            path = bfs_path(current_state[box], goal_position, current_state)
            for step in path[1:]:  # Skip the first step as it's the current position
                current_state[box] = step
                plan.append(current_state.copy())
    
    return plan

# Generate the plan
plan = generate_plan(initial_state, goal_state)
print(plan)