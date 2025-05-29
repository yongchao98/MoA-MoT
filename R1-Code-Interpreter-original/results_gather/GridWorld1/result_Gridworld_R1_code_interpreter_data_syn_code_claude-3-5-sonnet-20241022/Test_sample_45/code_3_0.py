import json
from collections import deque

# Initialize the problem
initial_pos = "C4,4"
goals = ['C4,2', 'C5,1', 'C5,4', 'C3,2', 'C3,3', 'C5,3']
obstacles = ['C1,2', 'C2,5', 'C5,2', 'C1,5', 'C3,5', 'C1,3', 'C1,4', 'C4,3']

# Adjacency dictionary
adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C2,5"], "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C1,5", "C3,5"], "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C2,5", "C4,5"], "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C3,5", "C5,5"], "C5,1": ["C5,2", "C4,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2"], "C5,3": ["C5,2", "C5,4", "C4,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4"], "C5,5": ["C5,4", "C4,5"]
}

def bfs_path(start, target, visited):
    if start == target:
        return [start]
    queue = deque([(start, [start])])
    seen = {start}
    
    while queue:
        current, path = queue.popleft()
        for next_pos in adjacency[current]:
            if next_pos not in seen and next_pos not in obstacles and next_pos not in visited:
                if next_pos == target:
                    return path + [next_pos]
                seen.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def find_complete_path():
    current = initial_pos
    full_path = [current]
    visited = set([current])
    remaining_goals = set(goals)
    
    # Remove initial position from goals if it's a goal
    if current in remaining_goals:
        remaining_goals.remove(current)
    
    while remaining_goals:
        found_path = None
        next_goal = None
        
        # Try each remaining goal
        for goal in remaining_goals:
            path = bfs_path(current, goal, set())
            if path and (found_path is None or len(path) < len(found_path)):
                found_path = path
                next_goal = goal
        
        if not found_path:
            return None
        
        # Add the path (excluding the start position)
        full_path.extend(found_path[1:])
        visited.update(found_path)
        current = next_goal
        remaining_goals.remove(next_goal)
    
    return full_path

# Find and print the path
path = find_complete_path()
print(f"<<<{json.dumps(path)}>>>")