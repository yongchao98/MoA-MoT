# Define the initial and goal states
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

# Define the adjacency list
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

# Function to find a path for a single box
def find_path(start, goal, occupied):
    from collections import deque
    
    queue = deque([(start, [start])])
    visited = set()
    
    while queue:
        current, path = queue.popleft()
        
        if current == goal:
            return path
        
        for neighbor in adjacency[current]:
            if neighbor not in visited and neighbor not in occupied:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    
    return []

# Generate the plan
def generate_plan(initial_state, goal_state):
    current_state = initial_state.copy()
    plan = [current_state.copy()]
    
    while current_state != goal_state:
        for box, goal in goal_state.items():
            if current_state[box] != goal:
                path = find_path(current_state[box], goal, set(current_state.values()))
                if len(path) > 1:
                    current_state[box] = path[1]
                    plan.append(current_state.copy())
                    break
    
    return plan

# Generate and print the plan
plan = generate_plan(initial_state, goal_state)
print(plan)