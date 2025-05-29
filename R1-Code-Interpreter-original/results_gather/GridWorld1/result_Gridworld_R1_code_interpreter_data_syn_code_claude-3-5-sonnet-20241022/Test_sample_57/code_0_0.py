from collections import deque
import json

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
        # Find nearest goal
        min_path = None
        nearest_goal = None
        
        for goal in remaining_goals:
            if current == goal:
                nearest_goal = goal
                min_path = [current]
                break
                
            temp_path = find_path(current, goal, adjacency, obstacles)
            if temp_path and (min_path is None or len(temp_path) < len(min_path)):
                min_path = temp_path
                nearest_goal = goal
        
        if min_path is None:
            return None
            
        # Add path to nearest goal (excluding the start position which is already in the path)
        path.extend(min_path[1:])
        remaining_goals.remove(nearest_goal)
        current = nearest_goal
    
    return path

# Define the problem
start = "C4,2"
goals = ["C1,3", "C2,4", "C3,4", "C3,2"]
obstacles = ["C3,5", "C4,3", "C3,1", "C1,5", "C1,4"]

# Define adjacency (use the provided adjacency dictionary)
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

# Find the complete path
result = find_complete_path(start, goals, adjacency, obstacles)

# Format and print the result
print(f"<<<{json.dumps(result)}>>>")