import json
from collections import deque

# Define the grid dimensions, obstacles, goals, and initial position
obstacles = {'C5,3', 'C1,2', 'C4,2', 'C2,7', 'C3,7', 'C2,2', 'C1,4', 'C6,2', 'C5,6', 'C4,5'}
goals = {'C1,1', 'C2,1', 'C5,5', 'C3,5', 'C6,4', 'C6,3', 'C5,1'}
initial_position = 'C3,3'

# Define adjacency list
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C2,7"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"],
    "C2,6": ["C2,5", "C2,7", "C1,6", "C3,6"],
    "C2,7": ["C2,6", "C1,7", "C3,7"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"],
    "C3,6": ["C3,5", "C3,7", "C2,6", "C4,6"],
    "C3,7": ["C3,6", "C2,7", "C4,7"],
    "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5", "C5,5"],
    "C4,6": ["C4,5", "C4,7", "C3,6", "C5,6"],
    "C4,7": ["C4,6", "C3,7", "C5,7"],
    "C5,1": ["C5,2", "C4,1", "C6,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2", "C6,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3", "C6,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4", "C6,4"],
    "C5,5": ["C5,4", "C5,6", "C4,5", "C6,5"],
    "C5,6": ["C5,5", "C5,7", "C4,6", "C6,6"],
    "C5,7": ["C5,6", "C4,7", "C6,7"],
    "C6,1": ["C6,2", "C5,1"],
    "C6,2": ["C6,1", "C6,3", "C5,2"],
    "C6,3": ["C6,2", "C6,4", "C5,3"],
    "C6,4": ["C6,3", "C6,5", "C5,4"],
    "C6,5": ["C6,4", "C6,6", "C5,5"],
    "C6,6": ["C6,5", "C6,7", "C5,6"],
    "C6,7": ["C6,6", "C5,7"]
}

# BFS to find a path visiting all goals
def bfs_find_path(start, goals, obstacles, adjacency):
    queue = deque([(start, [start], set())])  # (current_position, path, visited_goals)
    visited = set()

    while queue:
        current_position, path, visited_goals = queue.popleft()

        # If all goals are visited, return the path
        if visited_goals == goals:
            return path

        # Explore neighbors
        for neighbor in adjacency.get(current_position, []):
            if neighbor not in obstacles and (neighbor, frozenset(visited_goals)) not in visited:
                new_visited_goals = visited_goals.copy()
                if neighbor in goals:
                    new_visited_goals.add(neighbor)
                queue.append((neighbor, path + [neighbor], new_visited_goals))
                visited.add((neighbor, frozenset(new_visited_goals)))

    return None

# Find the path
path = bfs_find_path(initial_position, goals, obstacles, adjacency)

# Output the path as a JSON list
print(json.dumps(path))