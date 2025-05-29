import json
from collections import deque

def find_path(start, end, adjacency, obstacles):
    queue = deque([(start, [start])])
    visited = set([start])
    
    while queue:
        current, path = queue.popleft()
        if current == end:
            return path
            
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                visited.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def find_complete_path(start, goals, adjacency, obstacles):
    current = start
    remaining_goals = set(goals)
    path = [start]
    
    while remaining_goals:
        # Find the nearest unvisited goal
        shortest_path = None
        nearest_goal = None
        for goal in remaining_goals:
            if current == goal:
                temp_path = [current]
            else:
                temp_path = find_path(current, goal, adjacency, obstacles)
            if temp_path and (shortest_path is None or len(temp_path) < len(shortest_path)):
                shortest_path = temp_path
                nearest_goal = goal
        
        if shortest_path is None:
            return None  # No valid path found
            
        # Add the path to the nearest goal (excluding the start position if it's not the first path)
        path.extend(shortest_path[1:])
        current = nearest_goal
        remaining_goals.remove(nearest_goal)
    
    return path

# Define the problem
start = "C4,2"
goals = ["C1,3", "C2,2", "C4,1", "C2,3"]
obstacles = ["C4,3", "C3,1", "C3,3", "C2,1", "C1,2"]

# Define adjacency (using the provided adjacency dictionary)
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

# Find the complete path
path = find_complete_path(start, goals, adjacency, obstacles)

# Output the result
print(json.dumps(path))