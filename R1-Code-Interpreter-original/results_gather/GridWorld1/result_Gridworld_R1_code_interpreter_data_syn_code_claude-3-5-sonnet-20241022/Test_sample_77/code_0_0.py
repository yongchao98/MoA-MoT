from collections import deque
import json

def find_path(start, end, adjacency, visited_goals, obstacles):
    queue = deque([(start, [start])])
    seen = {start}
    
    while queue:
        current, path = queue.popleft()
        if current == end:
            return path
            
        for next_pos in adjacency[current]:
            if next_pos not in seen and next_pos not in obstacles:
                seen.add(next_pos)
                new_path = path + [next_pos]
                queue.append((next_pos, new_path))
    return None

def find_complete_path(start, goals, adjacency, obstacles):
    current_pos = start
    remaining_goals = set(goals)
    final_path = [start]
    visited_goals = set()
    
    while remaining_goals:
        # Find the nearest unvisited goal
        best_path = None
        best_goal = None
        best_length = float('inf')
        
        for goal in remaining_goals:
            path = find_path(current_pos, goal, adjacency, visited_goals, obstacles)
            if path and len(path) < best_length:
                best_length = len(path)
                best_path = path
                best_goal = goal
        
        if best_path is None:
            return None  # No valid path found
            
        # Add the path to the goal (excluding the start position)
        final_path.extend(best_path[1:])
        current_pos = best_goal
        remaining_goals.remove(best_goal)
        visited_goals.add(best_goal)
    
    return final_path

# Define the problem
start = "C2,4"
goals = ["C1,4", "C3,2", "C4,1", "C4,3"]
obstacles = ["C3,3", "C1,1", "C2,1", "C3,1", "C4,4"]

# Load adjacency data
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

# Find the solution
solution = find_complete_path(start, goals, adjacency, obstacles)
print(f"<<<{json.dumps(solution)}>>>")