import json
from collections import deque
from itertools import permutations

def find_path(adj, start, end, obstacles):
    queue = deque([(start, [start])])
    visited = set([start])
    
    while queue:
        current, path = queue.popleft()
        if current == end:
            return path
            
        for next_pos in adj[current]:
            if next_pos not in visited and next_pos not in obstacles:
                visited.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def find_complete_path(adj, start, goals, obstacles):
    # Find all possible paths between each pair of points
    paths = {}
    all_points = [start] + goals
    
    for p1 in all_points:
        paths[p1] = {}
        for p2 in all_points:
            if p1 != p2:
                path = find_path(adj, p1, p2, obstacles)
                if path:
                    paths[p1][p2] = path

    # Try all possible goal orderings
    best_path = None
    min_length = float('inf')
    
    for goal_order in permutations(goals):
        current = start
        full_path = [start]
        valid = True
        
        for goal in goal_order:
            if goal not in paths[current]:
                valid = False
                break
            path_segment = paths[current][goal][1:]  # Exclude start position
            full_path.extend(path_segment)
            current = goal
            
        if valid and len(full_path) < min_length:
            min_length = len(full_path)
            best_path = full_path
            
    return best_path

# Problem setup
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

start = "C3,2"
goals = ["C4,1", "C4,4", "C2,4", "C4,2"]
obstacles = ["C1,4", "C3,3", "C2,1", "C2,2", "C2,3"]

# Find the path
path = find_complete_path(adjacency, start, goals, obstacles)

# Output the result
print(f"<<<{json.dumps(path)}>>>")