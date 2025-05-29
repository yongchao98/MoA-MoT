import json
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

# Plan initialization
plan = [initial_state.copy()]

# Function to find a path for a single box using BFS
def bfs_path(start, goal, occupied):
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
current_state = initial_state.copy()
while current_state != goal_state:
    for box, current_position in current_state.items():
        if current_position != goal_state[box]:
            path = bfs_path(current_position, goal_state[box], set(current_state.values()))
            if path:
                for step in path[1:]:
                    current_state[box] = step
                    plan.append(current_state.copy())
                break

# Output the plan
print(json.dumps(plan, indent=2))