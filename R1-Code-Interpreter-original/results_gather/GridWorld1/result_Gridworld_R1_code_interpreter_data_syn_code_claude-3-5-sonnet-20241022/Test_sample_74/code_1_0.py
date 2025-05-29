import json
from collections import deque, defaultdict
import heapq

def find_shortest_path(start, target, adjacency, obstacles):
    if start in obstacles or target in obstacles:
        return None
    
    queue = [(0, start, [start])]
    visited = set()
    
    while queue:
        (cost, current, path) = heapq.heappop(queue)
        
        if current == target:
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
    
    while remaining_goals:
        # Find the nearest unvisited goal
        best_distance = float('inf')
        best_path = None
        best_goal = None
        
        for goal in remaining_goals:
            current_path = find_shortest_path(current, goal, adjacency, obstacles)
            if current_path is not None:
                distance = len(current_path)
                if distance < best_distance:
                    best_distance = distance
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
goals = ['C3,1', 'C2,1', 'C3,3', 'C1,1', 'C4,1', 'C5,2', 'C3,2']
obstacles = ['C5,3', 'C1,4', 'C1,2', 'C1,5', 'C5,4', 'C4,3', 'C2,2', 'C2,6', 'C4,6', 'C3,4']

# Adjacency dictionary
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

# Find the complete path
result = find_complete_path(start, goals, adjacency, obstacles)

# Format and print the result
if result:
    print(f"<<<{json.dumps(result)}>>>")
else:
    print("No valid path found")