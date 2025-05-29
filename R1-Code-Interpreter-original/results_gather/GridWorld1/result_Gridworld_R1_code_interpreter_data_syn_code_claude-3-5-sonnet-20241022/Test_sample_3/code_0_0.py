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
        min_dist = float('inf')
        next_goal = None
        next_path = None
        
        for goal in unvisited_goals:
            path = find_shortest_path(adjacency, current_pos, goal, obstacles)
            if path and len(path) < min_dist:
                min_dist = len(path)
                next_goal = goal
                next_path = path
        
        if next_goal is None:
            return None  # No valid path found
            
        # Add the path to this goal (excluding the start position)
        complete_path.extend(next_path[1:])
        unvisited_goals.remove(next_goal)
        current_pos = next_goal
    
    return complete_path

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

start = "C4,2"
goals = ['C5,1', 'C2,4', 'C3,4', 'C1,3', 'C2,2', 'C1,4']
obstacles = ['C3,3', 'C2,1', 'C5,2', 'C2,3', 'C5,3', 'C2,5', 'C1,1', 'C3,5']

# Find the complete path
path = find_complete_path(adjacency, start, goals, obstacles)

# Format the output
if path:
    print(f"<<<{json.dumps(path)}>>>")
else:
    print("No valid path found")