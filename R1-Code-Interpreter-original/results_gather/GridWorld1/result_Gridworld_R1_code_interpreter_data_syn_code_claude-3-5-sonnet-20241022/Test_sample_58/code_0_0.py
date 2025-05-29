import json
from collections import deque, defaultdict
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
    final_path = [start]
    
    while unvisited_goals:
        # Find the nearest unvisited goal
        best_path = None
        best_goal = None
        best_length = float('inf')
        
        for goal in unvisited_goals:
            path = find_shortest_path(adjacency, current_pos, goal, obstacles)
            if path and len(path) < best_length:
                best_length = len(path)
                best_path = path
                best_goal = goal
        
        if best_path is None:
            return None  # No valid path found
        
        # Add the path to the goal (excluding the start position)
        final_path.extend(best_path[1:])
        unvisited_goals.remove(best_goal)
        current_pos = best_goal
    
    return final_path

# Define the problem parameters
start = "C5,3"
goals = ['C5,4', 'C4,5', 'C1,3', 'C3,4', 'C3,3', 'C1,2']
obstacles = ['C1,4', 'C3,5', 'C4,3', 'C5,2', 'C5,1', 'C1,5', 'C2,1', 'C5,5']

# Load adjacency data
adjacency = {
    # ... (your adjacency dictionary)
}

# Find the complete path
path = find_complete_path(adjacency, start, goals, obstacles)

# Format and print the result
if path:
    result = json.dumps(path)
    print(f"<<<{result}>>>")
else:
    print("No valid path found")