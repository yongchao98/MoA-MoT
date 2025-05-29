from collections import deque

# Define the grid and parameters
initial_position = "C4,2"
goals = {"C2,4", "C3,3", "C4,1", "C1,2"}
obstacles = {"C4,3", "C3,4", "C4,4", "C1,1", "C2,3"}
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C2,4"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C1,4", "C3,4"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C2,4", "C4,4"],
    "C4,1": ["C4,2", "C3,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"],
    "C4,4": ["C4,3", "C3,4"]
}

# Remove obstacles from adjacency list
for obstacle in obstacles:
    if obstacle in adjacency:
        del adjacency[obstacle]
    for key in adjacency:
        if obstacle in adjacency[key]:
            adjacency[key].remove(obstacle)

# BFS to find a path that visits all goals
def bfs_path(start, goals):
    queue = deque([(start, [start], set())])
    while queue:
        position, path, visited_goals = queue.popleft()
        if visited_goals == goals:
            return path
        for neighbor in adjacency.get(position, []):
            if neighbor not in path:
                new_visited_goals = visited_goals | ({neighbor} if neighbor in goals else set())
                queue.append((neighbor, path + [neighbor], new_visited_goals))
    return []

# Find the path
path = bfs_path(initial_position, goals)
print(path)