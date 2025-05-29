import json
from collections import deque

def find_path(start, goals, obstacles, adjacency):
    def is_valid(cell):
        return cell not in obstacles
    
    def get_path_to_goal(start, target, visited):
        queue = deque([(start, [start])])
        seen = {start}
        
        while queue:
            current, path = queue.popleft()
            if current == target:
                return path
                
            for next_cell in adjacency[current]:
                if next_cell not in seen and is_valid(next_cell):
                    seen.add(next_cell)
                    queue.append((next_cell, path + [next_cell]))
        return None

    def find_complete_path():
        remaining_goals = set(goals)
        current = start
        final_path = [start]
        
        while remaining_goals:
            # Find closest remaining goal
            min_path = None
            closest_goal = None
            
            for goal in remaining_goals:
                path = get_path_to_goal(current, goal, set())
                if path is not None:
                    if min_path is None or len(path) < len(min_path):
                        min_path = path
                        closest_goal = goal
            
            if min_path is None:
                return None  # No valid path found
            
            # Add path to closest goal (excluding start position)
            final_path.extend(min_path[1:])
            remaining_goals.remove(closest_goal)
            current = closest_goal
            
        return final_path

    return find_complete_path()

# Set up the problem
start = "C4,4"
goals = ["C4,2", "C3,1", "C3,3", "C3,2"]
obstacles = ["C1,5", "C1,3", "C3,4", "C3,5", "C2,1"]

# Load adjacency list
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

# Find and print the path
path = find_path(start, goals, obstacles, adjacency)
print(f"<<<{json.dumps(path)}>>>")