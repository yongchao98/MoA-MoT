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

def find_complete_path(graph, start, goals, obstacles):
    unvisited_goals = set(goals)
    current = start
    complete_path = [start]
    
    while unvisited_goals:
        # Find the nearest unvisited goal
        min_path = None
        min_goal = None
        min_length = float('infinity')
        
        for goal in unvisited_goals:
            path = find_shortest_path(graph, current, goal, obstacles)
            if path and len(path) < min_length:
                min_length = len(path)
                min_path = path
                min_goal = goal
        
        if not min_path:
            return None
            
        # Add the path to the goal (excluding the start position)
        complete_path.extend(min_path[1:])
        unvisited_goals.remove(min_goal)
        current = min_goal
    
    return complete_path

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

start = "C2,2"
goals = ["C1,4", "C4,4", "C3,5", "C4,3", "C4,2", "C5,4"]
obstacles = ["C1,5", "C3,4", "C5,3", "C1,3", "C5,5", "C1,2", "C3,2", "C4,5"]

# Find the path
path = find_complete_path(graph, start, goals, obstacles)
print(f"<<<{json.dumps(path)}>>>")