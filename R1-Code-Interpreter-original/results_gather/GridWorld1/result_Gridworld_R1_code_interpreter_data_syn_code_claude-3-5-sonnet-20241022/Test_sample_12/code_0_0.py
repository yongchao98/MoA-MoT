import json
from collections import deque, defaultdict
from heapq import heappush, heappop

def heuristic(pos1, pos2):
    # Manhattan distance
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def find_path(start, end, adjacency, obstacles):
    if start == end:
        return [start]
    
    # A* search
    frontier = [(0, start, [start])]
    visited = set()
    
    while frontier:
        cost, current, path = heappop(frontier)
        
        if current == end:
            return path
            
        if current in visited:
            continue
            
        visited.add(current)
        
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                new_path = path + [next_pos]
                priority = len(new_path) + heuristic(next_pos, end)
                heappush(frontier, (priority, next_pos, new_path))
    
    return None

def find_complete_path():
    # Initialize data
    start = "C1,5"
    goals = ['C5,5', 'C1,6', 'C1,2', 'C2,2', 'C1,1', 'C5,2', 'C3,1']
    obstacles = ['C4,3', 'C4,5', 'C4,6', 'C2,6', 'C3,3', 'C4,2', 'C4,4', 'C2,5', 'C2,1', 'C6,1']
    
    # Load adjacency list
    adjacency = {
        # ... (the full adjacency dict as provided)
    }
    
    # Find nearest unvisited goal at each step
    current = start
    unvisited_goals = set(goals)
    final_path = [start]
    
    while unvisited_goals:
        # Find nearest goal
        min_dist = float('inf')
        next_goal = None
        best_path = None
        
        for goal in unvisited_goals:
            path = find_path(current, goal, adjacency, obstacles)
            if path and len(path) < min_dist:
                min_dist = len(path)
                next_goal = goal
                best_path = path
        
        if best_path is None:
            return None  # No valid path found
            
        # Add path to final path (excluding start position which is already included)
        final_path.extend(best_path[1:])
        current = next_goal
        unvisited_goals.remove(next_goal)
    
    return final_path

# Execute and format result
result = find_complete_path()
print(f"<<<{json.dumps(result)}>>>")