from collections import deque
import json
import heapq

def get_path(start, end, adjacency, obstacles):
    if start == end:
        return [start]
    
    # A* search
    frontier = [(0, start, [start])]
    visited = set()
    
    while frontier:
        cost, current, path = heapq.heappop(frontier)
        
        if current == end:
            return path
            
        if current in visited:
            continue
            
        visited.add(current)
        
        for next_pos in adjacency[current]:
            if next_pos not in obstacles and next_pos not in visited:
                new_path = path + [next_pos]
                # Simple Manhattan distance heuristic
                x1, y1 = map(int, next_pos.replace('C','').split(','))
                x2, y2 = map(int, end.replace('C','').split(','))
                h = abs(x1-x2) + abs(y1-y2)
                heapq.heappush(frontier, (len(new_path) + h, next_pos, new_path))
    
    return None

def find_complete_path(start, goals, adjacency, obstacles):
    current = start
    remaining_goals = set(goals)
    path = [start]
    
    while remaining_goals:
        # Find nearest unvisited goal
        best_goal = None
        best_path = None
        best_length = float('inf')
        
        for goal in remaining_goals:
            current_path = get_path(current, goal, adjacency, obstacles)
            if current_path and len(current_path) < best_length:
                best_length = len(current_path)
                best_path = current_path
                best_goal = goal
        
        if best_path is None:
            return None
        
        # Add path to goal (excluding start position)
        path.extend(best_path[1:])
        current = best_goal
        remaining_goals.remove(best_goal)
    
    return path

# Define the problem
start = "C2,1"
goals = ["C3,4", "C2,2", "C2,3", "C1,2"]
obstacles = ["C1,4", "C4,5", "C2,5", "C4,2", "C1,3"]

# Load adjacency data
adjacency = {
    # ... (your provided adjacency dict)
}

# Find solution
solution = find_complete_path(start, goals, adjacency, obstacles)

# Output solution
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")