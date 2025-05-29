from collections import deque
import json

def find_shortest_path(graph, start, end, obstacles):
    queue = deque([[start]])
    visited = set([start])
    
    while queue:
        path = queue.popleft()
        current = path[-1]
        
        if current == end:
            return path
            
        for neighbor in graph[current]:
            if neighbor not in visited and neighbor not in obstacles:
                queue.append(path + [neighbor])
                visited.add(neighbor)
    return None

def find_multi_goal_path(graph, start, goals, obstacles):
    current_pos = start
    remaining_goals = set(goals)
    final_path = [start]
    
    while remaining_goals:
        # Find the nearest unvisited goal
        shortest_dist = float('inf')
        next_goal = None
        next_path = None
        
        for goal in remaining_goals:
            path = find_shortest_path(graph, current_pos, goal, obstacles)
            if path and len(path) < shortest_dist:
                shortest_dist = len(path)
                next_goal = goal
                next_path = path
        
        if next_path is None:
            return None  # No valid path found
            
        # Add the path to the final path (excluding the start position)
        final_path.extend(next_path[1:])
        current_pos = next_goal
        remaining_goals.remove(next_goal)
    
    return final_path

# Define the problem
graph = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C2,4"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C1,4", "C3,4"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C2,4", "C4,4"],
    "C4,1": ["C4,2", "C3,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"],
    "C4,4": ["C4,3", "C3,4"]
}

start = "C2,4"
goals = ["C3,3", "C1,4", "C1,2", "C1,3"]
obstacles = ["C2,2", "C4,4", "C4,2", "C1,1", "C2,1"]

# Find the path
path = find_multi_goal_path(graph, start, goals, obstacles)

# Output the result
print(f"<<<{json.dumps(path)}>>>")