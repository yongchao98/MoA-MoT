import json
from collections import deque

def find_path(start, end, adjacency, visited_goals):
    queue = deque([(start, [start])])
    visited = set([start])
    
    while queue:
        current, path = queue.popleft()
        if current == end:
            return path
            
        for next_pos in adjacency[current]:
            if next_pos not in visited:
                visited.add(next_pos)
                new_path = path + [next_pos]
                queue.append((next_pos, new_path))
    return None

def find_complete_path(start, goals, adjacency, obstacles):
    remaining_goals = set(goals)
    current_pos = start
    final_path = [start]
    visited_goals = set()
    
    while remaining_goals:
        # Find shortest path to any remaining goal
        best_path = None
        best_goal = None
        best_length = float('inf')
        
        for goal in remaining_goals:
            path = find_path(current_pos, goal, adjacency, visited_goals)
            if path and len(path) < best_length:
                best_length = len(path)
                best_path = path
                best_goal = goal
        
        if not best_path:
            return None
            
        # Add the path (excluding the start position if it's not the first path)
        final_path.extend(best_path[1:])
        remaining_goals.remove(best_goal)
        visited_goals.add(best_goal)
        current_pos = best_goal
    
    return final_path

# Define the problem
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

start = "C2,2"
goals = ["C1,1", "C2,3", "C1,2", "C3,4"]
obstacles = ["C4,1", "C3,1", "C2,1", "C4,2", "C4,4"]

# Remove obstacle connections from adjacency
for obs in obstacles:
    if obs in adjacency:
        # Remove this position from its neighbors' adjacency lists
        for neighbor in adjacency[obs]:
            adjacency[neighbor].remove(obs)
        # Remove this position's adjacency list
        del adjacency[obs]

# Find the complete path
path = find_complete_path(start, goals, adjacency, obstacles)

# Output the result
print(f"<<<{json.dumps(path)}>>>")