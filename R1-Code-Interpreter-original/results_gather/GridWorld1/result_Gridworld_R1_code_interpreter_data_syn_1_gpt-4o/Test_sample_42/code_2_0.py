import json
from collections import deque

# Define the grid and constraints
obstacles = {'C2,3', 'C3,3', 'C3,2', 'C2,4', 'C3,4'}
goals = {'C1,2', 'C4,1', 'C1,3', 'C4,4'}
initial_position = 'C1,1'
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3"],
    "C1,3": ["C1,2", "C1,4"],
    "C1,4": ["C1,3"],
    "C2,1": ["C1,1", "C3,1"],
    "C2,2": ["C1,2"],
    "C3,1": ["C2,1", "C4,1"],
    "C4,1": ["C3,1", "C4,2"],
    "C4,2": ["C4,1", "C4,3"],
    "C4,3": ["C4,2", "C4,4"],
    "C4,4": ["C4,3"]
}

# BFS to find a path visiting all goals
def find_path():
    queue = deque([(initial_position, [initial_position], set())])
    while queue:
        current, path, visited_goals = queue.popleft()
        if visited_goals == goals:
            return path
        for neighbor in adjacency.get(current, []):
            if neighbor not in obstacles and neighbor not in path:
                new_visited_goals = visited_goals.copy()
                if neighbor in goals:
                    new_visited_goals.add(neighbor)
                queue.append((neighbor, path + [neighbor], new_visited_goals))
    return []

# Find the path
path = find_path()

# Output the path as a JSON list
print(json.dumps(path))