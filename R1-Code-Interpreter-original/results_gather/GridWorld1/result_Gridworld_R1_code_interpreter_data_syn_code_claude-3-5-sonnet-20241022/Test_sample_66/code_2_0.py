import json
from collections import deque
from itertools import permutations

def find_shortest_path(graph, start, end, obstacles):
    if start in obstacles or end in obstacles:
        return None
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        if current == end:
            return path
            
        for neighbor in graph[current]:
            if neighbor not in visited and neighbor not in obstacles:
                visited.add(neighbor)
                new_path = list(path)
                new_path.append(neighbor)
                queue.append((neighbor, new_path))
    return None

def merge_paths(paths):
    if not paths:
        return []
    result = list(paths[0])
    for path in paths[1:]:
        result.extend(path[1:])  # Skip the first element to avoid repetition
    return result

def find_optimal_path_with_ordering(graph, start, goals, obstacles):
    best_path = None
    min_length = float('infinity')
    
    # Try different orderings of goals
    for goal_order in permutations(goals):
        current_pos = start
        current_paths = []
        valid_sequence = True
        
        for goal in goal_order:
            path = find_shortest_path(graph, current_pos, goal, obstacles)
            if not path:
                valid_sequence = False
                break
            current_paths.append(path)
            current_pos = goal
        
        if valid_sequence:
            complete_path = merge_paths(current_paths)
            if len(complete_path) < min_length:
                min_length = len(complete_path)
                best_path = complete_path
                
    return best_path

# Input data
graph = {
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

start = "C1,2"
goals = ["C4,4", "C4,3", "C4,5", "C3,4", "C3,3", "C5,1"]
obstacles = ["C4,2", "C3,5", "C2,2", "C3,1", "C1,5", "C4,1", "C2,4", "C3,2"]

# Find the path
path = find_optimal_path_with_ordering(graph, start, goals, obstacles)

# Output the result
if path:
    print(f"<<<{json.dumps(path)}>>>")
else:
    print("No valid path found")