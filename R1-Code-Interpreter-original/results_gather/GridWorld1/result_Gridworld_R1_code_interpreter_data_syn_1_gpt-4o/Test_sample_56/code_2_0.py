import json

# Define the grid parameters
initial_position = "C3,4"
goals = {"C3,2", "C1,3", "C1,1", "C2,2"}
obstacles = {"C4,3", "C1,4", "C4,1", "C2,4", "C1,2"}
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

# DFS to find a path that visits all goals
def dfs(current, path, visited_goals):
    if visited_goals == goals:
        return path
    for neighbor in adjacency.get(current, []):
        if neighbor not in obstacles and neighbor not in path:
            new_visited_goals = visited_goals | ({neighbor} if neighbor in goals else set())
            result = dfs(neighbor, path + [neighbor], new_visited_goals)
            if result:
                return result
    return None

# Find the path
path = dfs(initial_position, [initial_position], set())

# Output the path as a JSON list
print(json.dumps(path))