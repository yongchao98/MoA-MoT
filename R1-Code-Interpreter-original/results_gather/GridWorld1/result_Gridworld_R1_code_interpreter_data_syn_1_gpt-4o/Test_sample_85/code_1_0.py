import json
from collections import deque

# Define the grid parameters
initial_position = "C4,2"
goals = {"C2,5", "C1,3", "C2,1", "C3,3"}
obstacles = {"C3,2", "C1,1", "C3,5", "C3,1", "C1,2"}
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
    "C4,1": ["C4,2", "C3,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C4,5": ["C4,4", "C3,5"]
}

# Function to find a path using BFS
def find_path():
    queue = deque([(initial_position, [initial_position], set())])
    
    while queue:
        position, current_path, visited_goals = queue.popleft()
        
        if position in goals:
            visited_goals = visited_goals | {position}
        
        if visited_goals == goals:
            return current_path
        
        for neighbor in adjacency.get(position, []):
            if neighbor not in obstacles and neighbor not in current_path:
                queue.append((neighbor, current_path + [neighbor], visited_goals))
    
    return []

# Find the path
path = find_path()

# Output the path as a JSON list
print(json.dumps(path))