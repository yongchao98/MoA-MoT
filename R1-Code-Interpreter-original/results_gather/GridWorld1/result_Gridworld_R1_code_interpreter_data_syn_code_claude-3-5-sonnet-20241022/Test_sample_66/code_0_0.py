import json
from collections import deque, defaultdict
import heapq

def find_shortest_path(graph, start, end, obstacles):
    if start in obstacles or end in obstacles:
        return None
    
    distances = {node: float('infinity') for node in graph}
    distances[start] = 0
    pq = [(0, start)]
    previous = {start: None}
    
    while pq:
        current_distance, current = heapq.heappop(pq)
        
        if current == end:
            path = []
            while current:
                path.append(current)
                current = previous[current]
            return path[::-1]
            
        if current_distance > distances[current]:
            continue
            
        for neighbor in graph[current]:
            if neighbor in obstacles:
                continue
                
            distance = current_distance + 1
            if distance < distances[neighbor]:
                distances[neighbor] = distance
                previous[neighbor] = current
                heapq.heappush(pq, (distance, neighbor))
    return None

def find_path_visiting_all_goals(graph, start, goals, obstacles):
    unvisited_goals = set(goals)
    current_pos = start
    final_path = [start]
    
    while unvisited_goals:
        # Find the nearest unvisited goal
        best_path = None
        best_goal = None
        min_distance = float('infinity')
        
        for goal in unvisited_goals:
            path = find_shortest_path(graph, current_pos, goal, obstacles)
            if path and len(path) < min_distance:
                min_distance = len(path)
                best_path = path
                best_goal = goal
        
        if not best_path:
            return None  # No valid path found
            
        # Add the path to the goal (excluding the start position)
        final_path.extend(best_path[1:])
        unvisited_goals.remove(best_goal)
        current_pos = best_goal
    
    return final_path

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
path = find_path_visiting_all_goals(graph, start, goals, obstacles)

# Output the result
if path:
    print(f"<<<{json.dumps(path)}>>>")
else:
    print("No valid path found")