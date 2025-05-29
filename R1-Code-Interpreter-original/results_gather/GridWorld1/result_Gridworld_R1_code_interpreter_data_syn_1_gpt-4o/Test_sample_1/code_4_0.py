import json
from collections import deque

# Define the grid and problem parameters
initial_position = "C2,5"
goals = {"C3,1", "C4,4", "C4,5", "C1,4"}
obstacles = {"C3,4", "C1,5", "C1,2", "C2,2", "C4,1"}
adjacency = {
    "C1,1": ["C2,1"],
    "C1,3": ["C1,4", "C2,3"],
    "C1,4": ["C1,3", "C2,4"],
    "C2,1": ["C1,1", "C3,1"],
    "C2,3": ["C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C1,4"],
    "C2,5": ["C2,4", "C3,5"],
    "C3,1": ["C3,2", "C2,1"],
    "C3,2": ["C3,1", "C3,3"],
    "C3,3": ["C3,2", "C4,3"],
    "C3,5": ["C2,5", "C4,5"],
    "C4,2": ["C4,3"],
    "C4,3": ["C4,2", "C4,4"],
    "C4,4": ["C4,3", "C4,5"],
    "C4,5": ["C4,4", "C3,5"]
}

# BFS to find a path visiting all goals
def find_path(initial, goals, adjacency):
    queue = deque([(initial, [initial], set())])  # (current_position, path, visited_goals)
    visited = set()  # To keep track of visited nodes to prevent cycles
    while queue:
        current, path, visited_goals = queue.popleft()
        if current in goals:
            visited_goals.add(current)
        if visited_goals == goals:
            return path
        for neighbor in adjacency.get(current, []):
            if neighbor not in path and neighbor not in obstacles and neighbor not in visited:  # Avoid cycles and obstacles
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor], visited_goals.copy()))
    return []

# Find the path
path = find_path(initial_position, goals, adjacency)
print(json.dumps(path))