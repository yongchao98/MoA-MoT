import json
from collections import deque
import heapq

def heuristic(a, b):
    # Manhattan distance
    x1, y1 = map(int, a.replace('C', '').split(','))
    x2, y2 = map(int, b.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def a_star(graph, start, goal, obstacles):
    frontier = [(0, start, [start])]
    visited = set()
    
    while frontier:
        f, current, path = heapq.heappop(frontier)
        
        if current == goal:
            return path
            
        if current in visited:
            continue
            
        visited.add(current)
        
        for next_pos in graph[current]:
            if next_pos not in visited and next_pos not in obstacles:
                new_path = path + [next_pos]
                priority = len(new_path) + heuristic(next_pos, goal)
                heapq.heappush(frontier, (priority, next_pos, new_path))
    
    return None

def find_path(graph, start, goals, obstacles):
    current = start
    path = [start]
    remaining_goals = set(goals)
    
    while remaining_goals:
        # Find nearest unvisited goal
        min_dist = float('inf')
        next_goal = None
        best_segment = None
        
        for goal in remaining_goals:
            segment = a_star(graph, current, goal, obstacles)
            if segment and len(segment) < min_dist:
                min_dist = len(segment)
                next_goal = goal
                best_segment = segment
        
        if not next_goal:
            return None
            
        # Add path segment (excluding start position which is already in path)
        path.extend(best_segment[1:])
        remaining_goals.remove(next_goal)
        current = next_goal
    
    return path

# Main execution
adjacency = {
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

start = "C2,3"
goals = ["C3,2", "C4,4", "C1,4", "C4,3", "C5,3", "C2,4"]
obstacles = ["C3,3", "C1,3", "C2,5", "C3,4", "C1,2", "C3,5", "C1,5", "C4,5"]

path = find_path(adjacency, start, goals, obstacles)
print(f"<<<{json.dumps(path)}>>>")