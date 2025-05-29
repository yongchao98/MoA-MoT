from collections import deque

# Define the grid and problem parameters
initial_position = "C2,5"
goals = {'C1,3', 'C4,1', 'C2,2', 'C4,4', 'C3,1', 'C1,4'}
obstacles = {'C3,2', 'C3,4', 'C3,5', 'C3,3', 'C5,4', 'C2,3', 'C5,1', 'C5,3'}
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C2,5"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C1,5", "C3,5"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C2,5", "C4,5"],
    "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C3,5", "C5,5"],
    "C5,1": ["C5,2", "C4,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4"],
    "C5,5": ["C5,4", "C4,5"]
}

# Remove obstacles from adjacency list
for obstacle in obstacles:
    if obstacle in adjacency:
        del adjacency[obstacle]
    for key in adjacency:
        if obstacle in adjacency[key]:
            adjacency[key].remove(obstacle)

# Breadth-First Search to find a path visiting all goals
def bfs_path(start, goals):
    queue = deque([(start, [start], set())])
    while queue:
        position, path, visited_goals = queue.popleft()
        if visited_goals == goals:
            return path
        for neighbor in adjacency.get(position, []):
            if neighbor not in path:
                new_visited_goals = visited_goals | {neighbor} if neighbor in goals else visited_goals
                queue.append((neighbor, path + [neighbor], new_visited_goals))
    return []

# Find the path
path = bfs_path(initial_position, goals)

# Output the path as a JSON list
import json
print(json.dumps(path))