import json
from collections import deque, defaultdict
import heapq

def find_shortest_path(start, end, adjacency, obstacles):
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

def find_complete_path(start, goals, adjacency, obstacles):
    current = start
    remaining_goals = set(goals)
    path = [start]
    
    if start in remaining_goals:
        remaining_goals.remove(start)
    
    while remaining_goals:
        # Find the nearest unvisited goal
        best_distance = float('inf')
        best_path = None
        best_goal = None
        
        for goal in remaining_goals:
            current_path = find_shortest_path(current, goal, adjacency, obstacles)
            if current_path is not None:
                if len(current_path) < best_distance:
                    best_distance = len(current_path)
                    best_path = current_path
                    best_goal = goal
        
        if best_path is None:
            return None  # No valid path found
            
        # Add the path to the goal (excluding the start position)
        path.extend(best_path[1:])
        current = best_goal
        remaining_goals.remove(best_goal)
    
    return path

# Problem setup
start = "C4,4"
goals = ['C2,4', 'C5,3', 'C3,4', 'C1,5', 'C6,6', 'C1,3', 'C4,2']
obstacles = ['C6,4', 'C3,5', 'C4,1', 'C2,1', 'C6,5', 'C1,4', 'C4,5', 'C5,1', 'C1,6', 'C1,2']

# Load adjacency data
adjacency = {
    # ... (your provided adjacency dictionary)
}

# Find the complete path
result = find_complete_path(start, goals, adjacency, obstacles)

# Format the output
if result:
    print(f"<<<{json.dumps(result)}>>>")
else:
    print("No valid path found")