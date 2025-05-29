import json
from collections import deque, defaultdict
import heapq

def find_shortest_path(adjacency, start, end, obstacles):
    if start in obstacles or end in obstacles:
        return None
    
    distances = {start: 0}
    previous = {start: None}
    pq = [(0, start)]
    visited = set()
    
    while pq:
        current_distance, current = heapq.heappop(pq)
        
        if current in visited:
            continue
            
        visited.add(current)
        
        if current == end:
            path = []
            while current:
                path.append(current)
                current = previous[current]
            return path[::-1]
            
        for neighbor in adjacency[current]:
            if neighbor in obstacles or neighbor in visited:
                continue
                
            distance = current_distance + 1
            if neighbor not in distances or distance < distances[neighbor]:
                distances[neighbor] = distance
                previous[neighbor] = current
                heapq.heappush(pq, (distance, neighbor))
    
    return None

def calculate_total_distance(path):
    return len(path) - 1

def find_complete_path(adjacency, start, goals, obstacles):
    current = start
    remaining_goals = set(goals)
    path = [start]
    visited_goals = set()
    
    while remaining_goals:
        # Find the nearest unvisited goal considering the overall path length
        min_total_distance = float('inf')
        next_goal = None
        next_path = None
        
        for goal in remaining_goals:
            current_path = find_shortest_path(adjacency, current, goal, obstacles)
            if current_path is not None:
                # Consider both distance to this goal and potential distances to other goals
                total_distance = calculate_total_distance(current_path)
                if total_distance < min_total_distance:
                    min_total_distance = total_distance
                    next_goal = goal
                    next_path = current_path
        
        if next_goal is None:
            return None  # No valid path found
            
        # Add the path to the next goal (excluding the start position)
        path.extend(next_path[1:])
        remaining_goals.remove(next_goal)
        visited_goals.add(next_goal)
        current = next_goal
    
    return path

# Input data
obstacles = ['C2,6', 'C5,2', 'C5,1', 'C5,4', 'C6,6', 'C4,1', 'C1,2', 'C1,1', 'C2,2', 'C5,7']
goals = ['C6,1', 'C2,3', 'C3,6', 'C1,4', 'C2,7', 'C6,3', 'C1,6']
start = 'C6,5'

# Define adjacency dictionary (same as before)
adjacency = {
    # ... (same adjacency dictionary as before)
}

# Find the complete path
result = find_complete_path(adjacency, start, goals, obstacles)

# Format the output
if result:
    output = json.dumps(result)
    print(f"<<<{output}>>>")
else:
    print("No valid path found")