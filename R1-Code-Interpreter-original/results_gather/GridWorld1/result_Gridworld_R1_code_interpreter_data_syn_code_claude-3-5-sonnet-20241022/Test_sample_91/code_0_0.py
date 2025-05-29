from collections import deque
import json

def find_path(start, target, adjacency, obstacles):
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        if current == target:
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
        chosen_goal = None
        
        for goal in remaining_goals:
            if current == goal:
                chosen_goal = goal
                min_path = [current]
                break
                
            temp_path = find_path(current, goal, adjacency, obstacles)
            if temp_path and (min_path is None or len(temp_path) < len(min_path)):
                min_path = temp_path
                chosen_goal = goal
        
        if min_path is None:
            return None
            
        # Add path to goal (excluding start position if it's not the first path)
        path.extend(min_path[1:] if len(path) > 0 else min_path)
        remaining_goals.remove(chosen_goal)
        current = chosen_goal
    
    return path

# Problem setup
start = "C4,2"
goals = ["C2,4", "C3,3", "C4,1", "C1,2"]
obstacles = ["C4,3", "C3,4", "C4,4", "C1,1", "C2,3"]

# Adjacency list
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

# Find solution
solution = find_complete_path(start, goals, adjacency, obstacles)
print(json.dumps(solution))