import json
from heapq import heappush, heappop
from itertools import permutations

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

def find_complete_path(start, goals, adjacency, obstacles):
    # First, find all valid paths between each pair of points
    points = [start] + goals
    paths = {}
    
    for p1 in points:
        for p2 in points:
            if p1 != p2:
                path = find_path(p1, p2, adjacency, obstacles)
                if path:
                    paths[(p1, p2)] = path

    # Try different goal orderings
    best_path = None
    min_length = float('inf')
    
    for goal_order in permutations(goals):
        current = start
        full_path = [start]
        valid = True
        
        for goal in goal_order:
            if (current, goal) not in paths:
                valid = False
                break
            path_segment = paths[(current, goal)][1:]  # Exclude start position
            full_path.extend(path_segment)
            current = goal
            
        if valid and len(full_path) < min_length:
            min_length = len(full_path)
            best_path = full_path
    
    return best_path

# Main execution
start = "C1,2"
goals = ["C1,4", "C4,4", "C2,1", "C3,1", "C3,4", "C3,5"]
obstacles = ["C4,2", "C2,2", "C3,2", "C5,2", "C4,3", "C3,3", "C5,3", "C2,3"]

# Load adjacency data
adjacency = {
    # ... (using the provided adjacency dictionary)
}

# Find the complete path
result = find_complete_path(start, goals, adjacency, obstacles)

# Output the result in the required format
print(f"<<<{json.dumps(result)}>>>")