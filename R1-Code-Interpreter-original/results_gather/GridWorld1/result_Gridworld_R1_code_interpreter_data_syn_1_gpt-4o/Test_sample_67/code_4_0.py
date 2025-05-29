import json

# Define the grid problem
initial_position = "C5,4"
goals = {'C4,5', 'C5,2', 'C1,5', 'C5,3', 'C4,2', 'C3,5', 'C3,6'}
obstacles = {'C3,4', 'C2,6', 'C4,3', 'C1,6', 'C2,5', 'C2,1', 'C2,4', 'C4,4', 'C4,1', 'C3,3'}
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
    "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5", "C5,5"],
    "C4,6": ["C4,5", "C3,6", "C5,6"],
    "C5,1": ["C5,2", "C4,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4"],
    "C5,5": ["C5,4", "C5,6", "C4,5"],
    "C5,6": ["C5,5", "C4,6"]
}

# DFS to find a path to a single goal
def dfs_path(initial, goal, obstacles, adjacency):
    stack = [(initial, [initial])]
    visited = set()
    while stack:
        current, path = stack.pop()
        if current == goal:
            return path
        if current not in visited:
            visited.add(current)
            for neighbor in adjacency.get(current, []):
                if neighbor not in obstacles and neighbor not in visited:
                    stack.append((neighbor, path + [neighbor]))
    return []

# Find a path visiting all goals
def find_full_path(initial, goals, obstacles, adjacency):
    current_position = initial
    full_path = [initial]
    remaining_goals = goals.copy()
    
    while remaining_goals:
        for goal in remaining_goals:
            path = dfs_path(current_position, goal, obstacles, adjacency)
            if path:
                full_path.extend(path[1:])  # Avoid duplicating the current position
                current_position = goal
                remaining_goals.remove(goal)
                break
    return full_path

# Find the full path
full_path = find_full_path(initial_position, goals, obstacles, adjacency)
print(json.dumps(full_path))