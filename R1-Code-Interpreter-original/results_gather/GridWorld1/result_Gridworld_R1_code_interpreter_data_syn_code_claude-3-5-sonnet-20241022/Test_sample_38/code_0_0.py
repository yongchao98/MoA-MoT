import json
from collections import deque

def find_shortest_path(start, end, adjacency, obstacles):
    if start in obstacles or end in obstacles:
        return None
    
    queue = deque([(start, [start])])
    visited = {start}
    
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
    unvisited_goals = set(goals)
    complete_path = [start]
    
    while unvisited_goals:
        # Find the nearest unvisited goal
        min_dist = float('inf')
        next_goal = None
        next_path = None
        
        for goal in unvisited_goals:
            path = find_shortest_path(current, goal, adjacency, obstacles)
            if path and len(path) < min_dist:
                min_dist = len(path)
                next_goal = goal
                next_path = path
        
        if next_goal is None:
            return None  # No valid path found
            
        # Add the path to this goal (excluding the start position)
        complete_path.extend(next_path[1:])
        unvisited_goals.remove(next_goal)
        current = next_goal
    
    return complete_path

# Problem setup
start = "C3,2"
goals = ['C2,1', 'C2,5', 'C1,5', 'C5,3', 'C5,1', 'C5,4']
obstacles = ['C5,2', 'C3,5', 'C2,4', 'C1,1', 'C5,5', 'C3,4', 'C3,1', 'C4,5']

# Load adjacency data
adjacency = {
    # ... (using the provided adjacency dictionary)
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
    "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C3,5", "C5,5"],
    "C5,1": ["C5,2", "C4,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4"],
    "C5,5": ["C5,4", "C4,5"]
}

# Find the complete path
path = find_complete_path(start, goals, adjacency, obstacles)

# Output the result
print(f"<<<{json.dumps(path)}>>>")