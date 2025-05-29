import json
from collections import deque

# Define the grid and problem parameters
initial_position = "C3,1"
goals = {'C5,4', 'C1,4', 'C2,4', 'C6,2', 'C2,6', 'C6,6', 'C3,5'}
obstacles = {'C1,3', 'C5,3', 'C3,3', 'C5,6', 'C5,2', 'C4,3', 'C1,2', 'C4,2', 'C6,3', 'C1,1'}
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
    "C5,1": ["C5,2", "C4,1", "C6,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2", "C6,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3", "C6,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4", "C6,4"],
    "C5,5": ["C5,4", "C5,6", "C4,5", "C6,5"],
    "C5,6": ["C5,5", "C4,6", "C6,6"],
    "C6,1": ["C6,2", "C5,1"],
    "C6,2": ["C6,1", "C6,3", "C5,2"],
    "C6,3": ["C6,2", "C6,4", "C5,3"],
    "C6,4": ["C6,3", "C6,5", "C5,4"],
    "C6,5": ["C6,4", "C6,6", "C5,5"],
    "C6,6": ["C6,5", "C5,6"]
}

# BFS to find a path to a single goal
def bfs_to_goal(start, goal, obstacles, adjacency):
    queue = deque([(start, [start])])
    visited = set()
    while queue:
        current, path = queue.popleft()
        if current == goal:
            return path
        for neighbor in adjacency.get(current, []):
            if neighbor not in obstacles and neighbor not in visited:
                queue.append((neighbor, path + [neighbor]))
                visited.add(neighbor)
    return []

# Find a path visiting all goals
def find_path(start, goals, obstacles, adjacency):
    current_position = start
    path = [start]
    remaining_goals = goals.copy()
    
    while remaining_goals:
        # Find the nearest goal
        nearest_goal = min(remaining_goals, key=lambda g: heuristic(current_position, {g}))
        # Find path to the nearest goal
        sub_path = bfs_to_goal(current_position, nearest_goal, obstacles, adjacency)
        if not sub_path:
            return []  # No path found
        # Update path and current position
        path.extend(sub_path[1:])  # Avoid duplicating the current position
        current_position = nearest_goal
        remaining_goals.remove(nearest_goal)
    
    return path

# Find the path
path = find_path(initial_position, goals, obstacles, adjacency)

# Output the path as a JSON list
print(json.dumps(path))