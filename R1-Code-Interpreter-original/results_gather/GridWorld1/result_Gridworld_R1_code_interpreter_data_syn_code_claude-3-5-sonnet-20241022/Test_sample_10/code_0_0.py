import json
from collections import deque
from itertools import permutations

def bfs_shortest_path(graph, start, end, obstacles):
    queue = deque([[start]])
    visited = {start}
    
    if start == end:
        return [start]
        
    while queue:
        path = queue.popleft()
        node = path[-1]
        
        for neighbor in graph[node]:
            if neighbor not in obstacles and neighbor not in visited:
                visited.add(neighbor)
                new_path = list(path)
                new_path.append(neighbor)
                
                if neighbor == end:
                    return new_path
                    
                queue.append(new_path)
    return None

def find_complete_path(graph, start, goals, obstacles):
    # First, get all valid goals (not blocked by obstacles)
    valid_goals = [g for g in goals if g not in obstacles]
    
    # Get all paths between each pair of points we need to visit
    paths = {}
    all_points = [start] + valid_goals
    
    for p1 in all_points:
        paths[p1] = {}
        for p2 in all_points:
            if p1 != p2:
                path = bfs_shortest_path(graph, p1, p2, set(obstacles))
                if path:
                    paths[p1][p2] = path

    # Try different goal orderings to find a valid complete path
    best_path = None
    min_length = float('inf')
    
    for perm in permutations(valid_goals):
        current_path = [start]
        current_pos = start
        valid_sequence = True
        
        for goal in perm:
            if goal not in paths[current_pos]:
                valid_sequence = False
                break
            # Add the path to the next goal (excluding the start point)
            current_path.extend(paths[current_pos][goal][1:])
            current_pos = goal
            
        if valid_sequence and len(current_path) < min_length:
            min_length = len(current_path)
            best_path = current_path
            
    return best_path

# Define the problem parameters
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

start = "C1,3"
goals = ["C5,1", "C4,3", "C5,4", "C1,4", "C3,4", "C3,2"]
obstacles = ["C5,2", "C5,3", "C4,5", "C5,5", "C3,5", "C2,4", "C1,5", "C2,5"]

# Find the solution path
solution = find_complete_path(adjacency, start, goals, obstacles)

# Output the solution in the required format
print(f"<<<{json.dumps(solution)}>>>")