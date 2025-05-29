import json
from collections import deque

# Define the grid dimensions, obstacles, goals, and initial position
grid_dimensions = (5, 6)
obstacles = {'C3,5', 'C1,5', 'C2,1', 'C2,3', 'C4,1', 'C5,5', 'C2,6', 'C4,3', 'C1,1', 'C1,6'}
goals = {'C5,1', 'C4,6', 'C2,2', 'C3,1', 'C3,4', 'C2,4', 'C5,3'}
initial_position = 'C4,2'

# Define adjacency list
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

# BFS to find a path that visits all goals
def bfs_find_path(start, goals, obstacles, adjacency):
    queue = deque([(start, [start], set())])  # (current_position, path, visited_goals)
    while queue:
        current_position, path, visited_goals = queue.popleft()
        if visited_goals == goals:
            return path
        for neighbor in adjacency.get(current_position, []):
            if neighbor not in obstacles and neighbor not in path:
                new_visited_goals = visited_goals | {neighbor} if neighbor in goals else visited_goals
                queue.append((neighbor, path + [neighbor], new_visited_goals))
    return None

# Find the path
path = bfs_find_path(initial_position, goals, obstacles, adjacency)

# Output the path as a JSON list
print(json.dumps(path))