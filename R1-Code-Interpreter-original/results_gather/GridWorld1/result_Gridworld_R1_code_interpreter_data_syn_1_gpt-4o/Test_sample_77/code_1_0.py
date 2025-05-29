from collections import deque

# Define the grid and problem parameters
initial_position = "C2,4"
goals = {"C1,4", "C3,2", "C4,1", "C4,3"}
obstacles = {"C3,3", "C1,1", "C2,1", "C3,1", "C4,4"}
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

# BFS to find a path visiting all goals
def find_path(initial, goals, obstacles, adjacency):
    queue = deque([(initial, [initial], set())])  # (current_position, path, visited_goals)
    visited = set()  # To keep track of visited nodes to prevent cycles
    while queue:
        current, path, visited_goals = queue.popleft()
        
        # Check if all goals are visited
        if visited_goals == goals:
            return path
        
        # Explore neighbors
        for neighbor in adjacency.get(current, []):
            if neighbor not in obstacles and neighbor not in visited:
                new_visited_goals = visited_goals | {neighbor} if neighbor in goals else visited_goals
                queue.append((neighbor, path + [neighbor], new_visited_goals))
                visited.add(neighbor)
    
    return []

# Find the path
path = find_path(initial_position, goals, obstacles, adjacency)
print(path)