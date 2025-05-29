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

# Function to find a path for a single box
def find_path(start, goal, adjacency):
    from collections import deque
    
    queue = deque([(start, [start])])
    visited = set()
    
    while queue:
        current, path = queue.popleft()
        
        if current == goal:
            return path
        
        for neighbor in adjacency[current]:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    
    return []

# Generate paths for each box
paths = {}
for box, start in initial_state.items():
    goal = goal_state[box]
    if start != goal:
        paths[box] = find_path(start, goal, adjacency)
    else:
        paths[box] = [start]

# Generate the plan
plan = [initial_state.copy()]
current_state = initial_state.copy()

# Iterate over the longest path length
max_steps = max(len(p) for p in paths.values())

for step in range(1, max_steps):
    new_state = current_state.copy()
    for box, path in paths.items():
        if step < len(path):
            new_state[box] = path[step]
    plan.append(new_state)
    current_state = new_state

# Ensure the final state matches the goal state
plan.append(goal_state)

# Print the plan
print(plan)