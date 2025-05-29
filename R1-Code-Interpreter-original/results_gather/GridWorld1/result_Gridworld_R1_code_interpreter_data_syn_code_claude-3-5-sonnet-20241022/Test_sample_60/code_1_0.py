import json
from collections import deque, defaultdict
import heapq

def find_shortest_path(graph, start, end, obstacles):
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
        
        for next_node in graph[current]:
            if next_node not in visited and next_node not in obstacles:
                heapq.heappush(queue, (cost + 1, next_node, path + [next_node]))
    
    return None

def find_complete_path(graph, start, goals, obstacles):
    current = start
    remaining_goals = set(goals)
    final_path = [start]
    visited_goals = set()
    
    while remaining_goals:
        # Find the nearest unvisited goal
        min_dist = float('inf')
        next_goal = None
        next_path = None
        
        for goal in remaining_goals:
            path = find_shortest_path(graph, current, goal, obstacles)
            if path and len(path) < min_dist:
                min_dist = len(path)
                next_goal = goal
                next_path = path
        
        if next_goal is None:
            return None  # No valid path found
            
        # Add the path to this goal (excluding the start position)
        final_path.extend(next_path[1:])
        current = next_goal
        remaining_goals.remove(next_goal)
        visited_goals.add(next_goal)
    
    return final_path

# Complete adjacency graph
graph = {
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

start = "C4,4"
goals = ['C3,3', 'C1,3', 'C3,4', 'C4,3', 'C2,3', 'C3,1', 'C1,6']
obstacles = ['C5,2', 'C1,5', 'C5,1', 'C1,2', 'C5,3', 'C4,6', 'C2,2', 'C3,5', 'C4,2', 'C1,4']

# Find the path
path = find_complete_path(graph, start, goals, obstacles)

# Output the result
if path:
    print(f"<<<{json.dumps(path)}>>>")
else:
    print("No valid path found")