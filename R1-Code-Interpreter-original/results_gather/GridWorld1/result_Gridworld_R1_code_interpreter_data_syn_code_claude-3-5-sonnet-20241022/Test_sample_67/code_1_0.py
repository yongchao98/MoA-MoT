import json
from heapq import heappush, heappop

def find_shortest_path(adjacency, start, end, obstacles):
    if start in obstacles or end in obstacles:
        return None
    
    queue = [(0, start, [start])]
    visited = set()
    
    while queue:
        (cost, current, path) = heappop(queue)
        
        if current == end:
            return path
            
        if current in visited:
            continue
            
        visited.add(current)
        
        for next_pos in adjacency.get(current, []):
            if next_pos not in visited and next_pos not in obstacles:
                new_path = path + [next_pos]
                heappush(queue, (cost + 1, next_pos, new_path))
    
    return None

def find_complete_path(adjacency, start, goals, obstacles):
    current_pos = start
    remaining_goals = set(goals)
    final_path = [start]
    
    while remaining_goals:
        best_path = None
        best_goal = None
        min_length = float('inf')
        
        for goal in remaining_goals:
            path = find_shortest_path(adjacency, current_pos, goal, obstacles)
            if path and len(path) < min_length:
                min_length = len(path)
                best_path = path
                best_goal = goal
        
        if not best_path:
            return None
        
        final_path.extend(best_path[1:])
        current_pos = best_goal
        remaining_goals.remove(best_goal)
    
    return final_path

# Complete adjacency dictionary
adjacency = {
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

start = "C5,4"
goals = ['C4,5', 'C5,2', 'C1,5', 'C5,3', 'C4,2', 'C3,5', 'C3,6']
obstacles = ['C3,4', 'C2,6', 'C4,3', 'C1,6', 'C2,5', 'C2,1', 'C2,4', 'C4,4', 'C4,1', 'C3,3']

# Find the path
path = find_complete_path(adjacency, start, goals, obstacles)

# Output the result
if path:
    print(f"<<<{json.dumps(path)}>>>")
else:
    print("No valid path found")