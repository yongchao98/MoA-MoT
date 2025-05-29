import json
from collections import deque, defaultdict
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

def find_path_through_goals(graph, start, goals, obstacles):
    current = start
    path = [start]
    remaining_goals = set(goals)
    
    while remaining_goals:
        min_dist = float('infinity')
        next_goal = None
        next_path = None
        
        for goal in remaining_goals:
            current_path = find_shortest_path(graph, current, goal, obstacles)
            if current_path is None:
                continue
            
            dist = len(current_path) - 1
            if dist < min_dist:
                min_dist = dist
                next_goal = goal
                next_path = current_path
        
        if next_goal is None:
            return None
            
        path.extend(next_path[1:])
        remaining_goals.remove(next_goal)
        current = next_goal
    
    return path

# Problem parameters
start = "C2,4"
goals = ['C4,3', 'C3,6', 'C3,3', 'C6,7', 'C4,4', 'C1,7', 'C2,5']
obstacles = ['C2,6', 'C3,4', 'C6,1', 'C1,2', 'C6,3', 'C5,6', 'C1,4', 'C5,1', 'C5,4', 'C4,6']

# Adjacency graph
adjacency = {
  "C1,1": ["C1,2", "C2,1"],
  "C1,2": ["C1,1", "C1,3", "C2,2"],
  "C1,3": ["C1,2", "C1,4", "C2,3"],
  "C1,4": ["C1,3", "C1,5", "C2,4"],
  "C1,5": ["C1,4", "C1,6", "C2,5"],
  "C1,6": ["C1,5", "C1,7", "C2,6"],
  "C1,7": ["C1,6", "C2,7"],
  "C2,1": ["C2,2", "C1,1", "C3,1"],
  "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
  "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
  "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
  "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"],
  "C2,6": ["C2,5", "C2,7", "C1,6", "C3,6"],
  "C2,7": ["C2,6", "C1,7", "C3,7"],
  "C3,1": ["C3,2", "C2,1", "C4,1"],
  "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
  "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
  "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
  "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"],
  "C3,6": ["C3,5", "C3,7", "C2,6", "C4,6"],
  "C3,7": ["C3,6", "C2,7", "C4,7"],
  "C4,1": ["C4,2", "C3,1", "C5,1"],
  "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
  "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
  "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
  "C4,5": ["C4,4", "C4,6", "C3,5", "C5,5"],
  "C4,6": ["C4,5", "C4,7", "C3,6", "C5,6"],
  "C4,7": ["C4,6", "C3,7", "C5,7"],
  "C5,1": ["C5,2", "C4,1", "C6,1"],
  "C5,2": ["C5,1", "C5,3", "C4,2", "C6,2"],
  "C5,3": ["C5,2", "C5,4", "C4,3", "C6,3"],
  "C5,4": ["C5,3", "C5,5", "C4,4", "C6,4"],
  "C5,5": ["C5,4", "C5,6", "C4,5", "C6,5"],
  "C5,6": ["C5,5", "C5,7", "C4,6", "C6,6"],
  "C5,7": ["C5,6", "C4,7", "C6,7"],
  "C6,1": ["C6,2", "C5,1"],
  "C6,2": ["C6,1", "C6,3", "C5,2"],
  "C6,3": ["C6,2", "C6,4", "C5,3"],
  "C6,4": ["C6,3", "C6,5", "C5,4"],
  "C6,5": ["C6,4", "C6,6", "C5,5"],
  "C6,6": ["C6,5", "C6,7", "C5,6"],
  "C6,7": ["C6,6", "C5,7"]
}

# Find the solution path
solution_path = find_path_through_goals(adjacency, start, goals, obstacles)

# Format and print the result
if solution_path:
    result = json.dumps(solution_path)
    print(f"<<<{result}>>>")
else:
    print("No valid path found")