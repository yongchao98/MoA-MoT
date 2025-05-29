import json
from heapq import heappush, heappop
from itertools import permutations

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def a_star(start, goal, adjacency, obstacles):
    if start in obstacles or goal in obstacles:
        return None
    
    frontier = [(0, start, [start])]
    visited = set()
    
    while frontier:
        f_score, current, path = heappop(frontier)
        
        if current == goal:
            return path
            
        if current in visited:
            continue
            
        visited.add(current)
        
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                new_path = path + [next_pos]
                h_score = manhattan_distance(next_pos, goal)
                g_score = len(new_path) - 1
                f_score = g_score + h_score
                heappush(frontier, (f_score, next_pos, new_path))
    
    return None

def find_complete_path(start, goals, adjacency, obstacles):
    # Find paths between all pairs of points we need to visit
    points = [start] + goals
    paths = {}
    
    for p1 in points:
        for p2 in points:
            if p1 != p2:
                path = a_star(p1, p2, adjacency, obstacles)
                if path:
                    paths[(p1, p2)] = path

    # Try different goal orderings
    best_path = None
    min_length = float('inf')
    
    for goal_order in permutations(goals):
        current = start
        total_path = [start]
        valid_sequence = True
        
        for goal in goal_order:
            if (current, goal) not in paths:
                valid_sequence = False
                break
            path_segment = paths[(current, goal)][1:]  # Exclude start position
            total_path.extend(path_segment)
            current = goal
            
        if valid_sequence and len(total_path) < min_length:
            min_length = len(total_path)
            best_path = total_path
    
    return best_path

# Main execution
start = "C4,2"
goals = ['C5,1', 'C2,4', 'C3,4', 'C1,3', 'C2,2', 'C1,4']
obstacles = ['C3,3', 'C2,1', 'C5,2', 'C2,3', 'C5,3', 'C2,5', 'C1,1', 'C3,5']

# Load adjacency data
adjacency = {
    # ... (your provided adjacency dict)
}

# Find the path
path = find_complete_path(start, goals, adjacency, obstacles)

# Output the result
print(f"<<<{json.dumps(path)}>>>")