import json
from collections import deque
from itertools import permutations

def bfs_shortest_path(adj, start, end, obstacles):
    if start == end:
        return [start]
    if start in obstacles or end in obstacles:
        return None
    
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        current = path[-1]
        
        for next_pos in adj[current]:
            if next_pos not in visited and next_pos not in obstacles:
                if next_pos == end:
                    return path + [next_pos]
                visited.add(next_pos)
                queue.append(path + [next_pos])
    return None

def find_complete_path(adj, start, goals, obstacles):
    # Create distance matrix between all points (including start and goals)
    points = [start] + goals
    n = len(points)
    distances = {}
    paths = {}
    
    for i in range(n):
        for j in range(n):
            if i != j:
                path = bfs_shortest_path(adj, points[i], points[j], obstacles)
                if path:
                    distances[(points[i], points[j])] = len(path) - 1
                    paths[(points[i], points[j])] = path
                else:
                    distances[(points[i], points[j])] = float('inf')
                    paths[(points[i], points[j])] = None
    
    # Try all possible goal permutations
    best_length = float('inf')
    best_path = None
    
    for perm in permutations(goals):
        current_length = 0
        valid_path = True
        current_path = [start]
        current_pos = start
        
        for next_goal in perm:
            if distances[(current_pos, next_goal)] == float('inf'):
                valid_path = False
                break
            current_length += distances[(current_pos, next_goal)]
            if len(current_path) > 1:
                current_path.pop()  # Remove duplicate position
            current_path.extend(paths[(current_pos, next_goal)])
            current_pos = next_goal
            
        if valid_path and current_length < best_length:
            best_length = current_length
            best_path = current_path
    
    return best_path

# Problem setup
adj = {
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

start = "C2,5"
goals = ["C1,3", "C4,1", "C2,2", "C4,4", "C3,1", "C1,4"]
obstacles = ["C3,2", "C3,4", "C3,5", "C3,3", "C5,4", "C2,3", "C5,1", "C5,3"]

# Find the solution
path = find_complete_path(adj, start, goals, obstacles)
print(f"<<<{json.dumps(path)}>>>")