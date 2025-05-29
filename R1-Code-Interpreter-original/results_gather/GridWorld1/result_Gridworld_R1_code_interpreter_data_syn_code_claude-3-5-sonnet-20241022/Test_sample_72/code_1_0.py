import json
from collections import deque
import heapq

def find_shortest_path(graph, start, end, obstacles):
    if start in obstacles or end in obstacles:
        return None
    
    distances = {cell: float('infinity') for cell in graph}
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
    current = start
    remaining_goals = set(goals)
    path = [start]
    
    while remaining_goals:
        # Find the nearest unvisited goal
        best_distance = float('infinity')
        best_goal = None
        best_path = None
        
        for goal in remaining_goals:
            current_path = find_shortest_path(graph, current, goal, obstacles)
            if current_path is not None:
                distance = len(current_path) - 1
                if distance < best_distance:
                    best_distance = distance
                    best_goal = goal
                    best_path = current_path
        
        if best_goal is None:
            return None  # No valid path found
            
        # Add the path to this goal (excluding the start position)
        path.extend(best_path[1:])
        current = best_goal
        remaining_goals.remove(best_goal)
    
    return path

# Input data with complete adjacency graph
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

start = "C5,5"
goals = ['C3,3', 'C5,4', 'C4,1', 'C4,4', 'C1,5', 'C2,2']
obstacles = ['C5,3', 'C1,3', 'C3,4', 'C1,1', 'C2,4', 'C4,2', 'C1,2', 'C5,1']

# Find the path
path = find_complete_path(graph, start, goals, obstacles)

# Output the result
if path:
    print(f"<<<{json.dumps(path)}>>>")
else:
    print("No valid path found")