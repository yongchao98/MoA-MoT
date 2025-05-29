import json
from collections import deque
import heapq

def find_shortest_path(adjacency, start, end, obstacles):
    if start in obstacles or end in obstacles:
        return None
    
    queue = [(0, start, [start])]
    visited = set()
    
    while queue:
        (cost, current, path) = heapq.heappop(queue)
        
        if current == end:
            return path
            
        if current in visited:
            continue
            
        visited.add(current)
        
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                new_path = path + [next_pos]
                heapq.heappush(queue, (cost + 1, next_pos, new_path))
    
    return None

def find_complete_path(adjacency, start, goals, obstacles):
    unvisited_goals = set(goals)
    current_pos = start
    complete_path = [start]
    
    while unvisited_goals:
        # Find the nearest unvisited goal
        best_path = None
        best_goal = None
        min_length = float('inf')
        
        for goal in unvisited_goals:
            path = find_shortest_path(adjacency, current_pos, goal, obstacles)
            if path and len(path) < min_length:
                min_length = len(path)
                best_path = path
                best_goal = goal
        
        if not best_path:
            return None  # No valid path found
            
        # Add the path to the goal (excluding the start position)
        complete_path.extend(best_path[1:])
        unvisited_goals.remove(best_goal)
        current_pos = best_goal
    
    return complete_path

# Input data
adjacency = {
    # ... (your provided adjacency dict)
}

start = "C2,6"
goals = ['C2,1', 'C2,3', 'C3,4', 'C5,2', 'C5,6', 'C4,3', 'C6,4']
obstacles = ['C1,3', 'C5,5', 'C2,5', 'C5,3', 'C3,3', 'C6,5', 'C1,6', 'C1,2', 'C3,5', 'C6,2']

# Find the path
path = find_complete_path(adjacency, start, goals, obstacles)

# Output the result
if path:
    print(f"<<<{json.dumps(path)}>>>")
else:
    print("No valid path found")