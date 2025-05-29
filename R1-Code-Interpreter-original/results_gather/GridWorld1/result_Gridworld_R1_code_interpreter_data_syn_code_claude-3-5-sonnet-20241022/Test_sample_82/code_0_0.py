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
    complete_path = [start]
    
    while unvisited_goals:
        # Find the nearest unvisited goal
        min_distance = float('inf')
        next_goal = None
        next_path = None
        
        for goal in unvisited_goals:
            path = find_shortest_path(adjacency, current_pos, goal, obstacles)
            if path and len(path) < min_distance:
                min_distance = len(path)
                next_goal = goal
                next_path = path
        
        if next_path is None:
            return None  # No valid path found
        
        # Add the path to the goal (excluding the start position)
        complete_path.extend(next_path[1:])
        unvisited_goals.remove(next_goal)
        current_pos = next_goal
    
    return complete_path

# Define the problem parameters
start = "C3,1"
goals = ['C5,4', 'C1,4', 'C2,4', 'C6,2', 'C2,6', 'C6,6', 'C3,5']
obstacles = ['C1,3', 'C5,3', 'C3,3', 'C4,3', 'C5,6', 'C5,2', 'C1,2', 'C4,2', 'C6,3', 'C1,1']

# Load adjacency data
adjacency = {
    # ... (using the provided adjacency dictionary)
}

# Find the complete path
path = find_complete_path(adjacency, start, goals, obstacles)

# Format the output
if path:
    output = json.dumps(path)
    print(f"<<<{output}>>>")
else:
    print("No valid path found")