from collections import deque
import json

# Setup the problem
initial_pos = "C1,5"
goals = ['C2,2', 'C1,4', 'C4,1', 'C4,3', 'C5,4', 'C2,5', 'C1,1']
obstacles = ['C3,4', 'C5,2', 'C5,5', 'C3,6', 'C4,2', 'C2,1', 'C1,6', 'C4,4', 'C2,6', 'C5,1']

# The adjacency dictionary (as provided)
adjacency = {
    # ... (using the provided adjacency dictionary)
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"],
    "C2,6": ["C2,5", "C1,6", "C3,6"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"],
    "C3,6": ["C3,5", "C2,6", "C4,6"],
    "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5", "C5,5"],
    "C4,6": ["C4,5", "C3,6", "C5,6"],
    "C5,1": ["C5,2", "C4,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4"],
    "C5,5": ["C5,4", "C5,6", "C4,5"],
    "C5,6": ["C5,5", "C4,6"]
}

def find_path(start, end, obstacles):
    if start == end:
        return [start]
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                if next_pos == end:
                    return path + [next_pos]
                visited.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def find_complete_path():
    current_pos = initial_pos
    unvisited_goals = goals.copy()
    complete_path = [current_pos]
    
    while unvisited_goals:
        # Find the closest unvisited goal
        min_dist = float('inf')
        next_goal = None
        best_path = None
        
        for goal in unvisited_goals:
            path = find_path(current_pos, goal, obstacles)
            if path and len(path) < min_dist:
                min_dist = len(path)
                next_goal = goal
                best_path = path
        
        if not best_path:
            return None  # No valid path found
            
        # Add the path to the complete path (excluding the first position as it's already included)
        complete_path.extend(best_path[1:])
        current_pos = next_goal
        unvisited_goals.remove(next_goal)
    
    return complete_path

# Find and print the solution
solution = find_complete_path()
print(f"<<<{json.dumps(solution)}>>>")