from collections import deque

# Define the grid and problem parameters
initial_position = "C2,2"
goals = set(['C1,4', 'C4,4', 'C3,5', 'C4,3', 'C4,2', 'C5,4'])
obstacles = set(['C1,5', 'C3,4', 'C5,3', 'C1,3', 'C5,5', 'C1,2', 'C3,2', 'C4,5'])

# Define adjacency list
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

# Filter adjacency list to exclude obstacles
filtered_adjacency = {k: [n for n in v if n not in obstacles] for k, v in adjacency.items()}

# BFS to find the path
def bfs_path(start, goals, adjacency):
    queue = deque([(start, [start], set())])  # (current_position, path, visited_goals)
    visited = set()  # To keep track of visited nodes to prevent cycles
    while queue:
        current, path, visited_goals = queue.popleft()
        if visited_goals == goals:
            return path
        if current in visited:
            continue
        visited.add(current)
        for neighbor in adjacency.get(current, []):
            if neighbor not in path:
                new_visited_goals = visited_goals.copy()
                if neighbor in goals:
                    new_visited_goals.add(neighbor)
                queue.append((neighbor, path + [neighbor], new_visited_goals))
    return []

# Find the path
path = bfs_path(initial_position, goals, filtered_adjacency)
print(path)