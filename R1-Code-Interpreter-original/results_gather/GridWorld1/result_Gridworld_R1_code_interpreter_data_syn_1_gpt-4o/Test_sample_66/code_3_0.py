import json
from collections import deque

# Define the grid dimensions, obstacles, goals, and initial position
obstacles = {'C4,2', 'C3,5', 'C2,2', 'C3,1', 'C1,5', 'C4,1', 'C2,4', 'C3,2'}
goals = {'C4,4', 'C4,3', 'C4,5', 'C3,4', 'C3,3', 'C5,1'}
initial_position = 'C1,2'

# Define the adjacency list, excluding obstacles
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3"],
    "C2,1": ["C1,1"],
    "C2,3": ["C1,3", "C3,3"],
    "C3,3": ["C2,3", "C4,3"],
    "C3,4": ["C4,4"],
    "C4,3": ["C3,3", "C4,4"],
    "C4,4": ["C3,4", "C4,5"],
    "C4,5": ["C4,4"],
    "C5,1": [],
}

# BFS to find a path that visits all goals
def bfs_path(start, goals, adjacency):
    queue = deque([(start, [start], set())])
    while queue:
        current, path, visited_goals = queue.popleft()
        if visited_goals == goals:
            return path
        for neighbor in adjacency.get(current, []):
            if neighbor not in path:
                new_visited_goals = visited_goals | {neighbor} if neighbor in goals else visited_goals
                queue.append((neighbor, path + [neighbor], new_visited_goals))
    return []

# Find the path
path = bfs_path(initial_position, goals, adjacency)

# Output the path as a JSON list
print(json.dumps(path))