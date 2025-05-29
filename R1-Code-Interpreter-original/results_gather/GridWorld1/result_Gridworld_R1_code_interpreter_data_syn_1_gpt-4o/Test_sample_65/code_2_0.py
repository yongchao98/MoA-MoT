import json

# Define the grid and obstacles
obstacles = {'C4,2', 'C3,5', 'C2,2', 'C3,1', 'C1,5', 'C4,1', 'C2,4', 'C3,2'}
goals = {'C4,4', 'C4,3', 'C4,5', 'C3,4', 'C3,3', 'C5,1'}
initial_position = 'C1,2'

# Define adjacency list
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
    "C5,1": []
}

# DFS to find a path visiting all goals
def dfs_path(position, path, visited_goals):
    if visited_goals == goals:
        return path
    for neighbor in adjacency.get(position, []):
        if neighbor not in path:
            new_visited_goals = visited_goals | {neighbor} if neighbor in goals else visited_goals
            result = dfs_path(neighbor, path + [neighbor], new_visited_goals)
            if result:
                return result
    return None

# Find the path
path = dfs_path(initial_position, [initial_position], set())

# Output the path as JSON
print(json.dumps(path))