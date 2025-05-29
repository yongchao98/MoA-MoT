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
        # Find the nearest unvisited goal
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
            return None  # No valid path found
            
        # Add the path to the next goal (excluding the start position)
        path.extend(next_path[1:])
        remaining_goals.remove(next_goal)
        current = next_goal
    
    return path

# Define the problem parameters
start = "C2,4"
goals = ['C4,3', 'C3,6', 'C3,3', 'C6,7', 'C4,4', 'C1,7', 'C2,5']
obstacles = ['C2,6', 'C3,4', 'C6,1', 'C1,2', 'C6,3', 'C5,6', 'C1,4', 'C5,1', 'C5,4', 'C4,6']

# Load the adjacency graph
adjacency = {
    # ... (the adjacency dict as provided in the problem)
}

# Find the solution path
solution_path = find_path_through_goals(adjacency, start, goals, obstacles)

# Format and print the result
if solution_path:
    result = json.dumps(solution_path)
    print(f"<<<{result}>>>")
else:
    print("No valid path found")