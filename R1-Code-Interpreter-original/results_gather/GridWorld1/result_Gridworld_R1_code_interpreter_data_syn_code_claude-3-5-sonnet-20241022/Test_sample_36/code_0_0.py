from collections import deque
import json

def get_manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def find_path(start, end, adjacency, obstacles):
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

def find_complete_path(start, goals, adjacency, obstacles):
    current_pos = start
    remaining_goals = goals.copy()
    complete_path = [start]
    
    while remaining_goals:
        # Find the nearest unvisited goal
        min_dist = float('inf')
        next_goal = None
        next_path = None
        
        for goal in remaining_goals:
            path = find_path(current_pos, goal, adjacency, obstacles)
            if path:
                dist = len(path) - 1
                if dist < min_dist:
                    min_dist = dist
                    next_goal = goal
                    next_path = path
        
        if next_goal is None:
            return None  # No valid path found
            
        # Add the path to this goal (excluding the start position)
        complete_path.extend(next_path[1:])
        current_pos = next_goal
        remaining_goals.remove(next_goal)
    
    return complete_path

# Problem setup
start = "C4,2"
goals = ["C2,1", "C1,4", "C3,1", "C1,2"]
obstacles = ["C1,1", "C3,4", "C2,4", "C1,5", "C3,5"]

# Adjacency dictionary (provided in the problem)
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

# Find the solution
solution = find_complete_path(start, goals, adjacency, obstacles)
print(f"<<<{json.dumps(solution)}>>>")