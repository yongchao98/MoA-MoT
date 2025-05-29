from collections import deque

# Define the problem parameters
initial_position = "C2,4"
goals = {"C3,3", "C1,4", "C1,2", "C1,3"}
obstacles = {"C2,2", "C4,4", "C4,2", "C1,1", "C2,1"}
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

# Remove obstacles from adjacency
for obstacle in obstacles:
    if obstacle in adjacency:
        del adjacency[obstacle]

# Remove obstacle neighbors
for key in adjacency:
    adjacency[key] = [neighbor for neighbor in adjacency[key] if neighbor not in obstacles]

# BFS to find a path visiting all goals
def find_path(initial, goals, adjacency):
    queue = deque([(initial, [initial], set())])  # (current_position, path, visited_goals)
    while queue:
        current, path, visited_goals = queue.popleft()
        if visited_goals == goals:
            return path
        for neighbor in adjacency.get(current, []):
            if neighbor not in path:  # Avoid cycles
                new_visited_goals = visited_goals | ({neighbor} if neighbor in goals else set())
                queue.append((neighbor, path + [neighbor], new_visited_goals))
    return []

# Find the path
path = find_path(initial_position, goals, adjacency)
print(path)