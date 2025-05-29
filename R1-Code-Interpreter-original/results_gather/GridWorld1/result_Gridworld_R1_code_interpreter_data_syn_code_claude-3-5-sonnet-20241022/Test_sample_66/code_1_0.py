import json
from collections import deque

def find_shortest_path(graph, start, end, obstacles):
    if start in obstacles or end in obstacles:
        return None
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        if current == end:
            return path
            
        for neighbor in graph[current]:
            if neighbor not in visited and neighbor not in obstacles:
                visited.add(neighbor)
                new_path = list(path)
                new_path.append(neighbor)
                queue.append((neighbor, new_path))
    return None

def find_optimal_path(graph, start, goals, obstacles):
    current = start
    path = [start]
    unvisited_goals = set(goals)
    visited_positions = {start}
    
    while unvisited_goals:
        best_path = None
        best_goal = None
        min_length = float('infinity')
        
        for goal in unvisited_goals:
            temp_path = find_shortest_path(graph, current, goal, obstacles)
            if temp_path and len(temp_path) < min_length:
                min_length = len(temp_path)
                best_path = temp_path
                best_goal = goal
        
        if not best_path:
            return None
        
        # Add the path to the goal (excluding the start position)
        for pos in best_path[1:]:
            path.append(pos)
            visited_positions.add(pos)
            # Mark any goals we pass through as visited
            if pos in unvisited_goals:
                unvisited_goals.remove(pos)
        
        current = best_goal
        if best_goal in unvisited_goals:
            unvisited_goals.remove(best_goal)
    
    return path

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
path = find_optimal_path(graph, start, goals, obstacles)

# Output the result
if path:
    print(f"<<<{json.dumps(path)}>>>")
else:
    print("No valid path found")