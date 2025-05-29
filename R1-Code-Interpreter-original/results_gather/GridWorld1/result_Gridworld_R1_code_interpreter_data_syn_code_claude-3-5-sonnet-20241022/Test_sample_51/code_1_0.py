import json
from heapq import heappush, heappop
from itertools import permutations

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def a_star_search(start, goal, adjacency, obstacles):
    if start in obstacles or goal in obstacles:
        return None
    
    frontier = [(0, start, [start])]
    visited = set()
    
    while frontier:
        cost, current, path = heappop(frontier)
        
        if current == goal:
            return path
            
        if current in visited:
            continue
            
        visited.add(current)
        
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                new_path = path + [next_pos]
                priority = len(new_path) + manhattan_distance(next_pos, goal)
                heappush(frontier, (priority, next_pos, new_path))
    
    return None

def find_complete_path(start, goals, adjacency, obstacles):
    # Find paths between all pairs of points (start and goals)
    all_points = [start] + goals
    paths = {}
    
    for p1 in all_points:
        for p2 in all_points:
            if p1 != p2:
                path = a_star_search(p1, p2, adjacency, obstacles)
                if path:
                    paths[(p1, p2)] = path

    # Try different goal orderings
    best_path = None
    min_length = float('inf')
    
    # Try all permutations of goals
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
start = "C2,5"
goals = ['C1,4', 'C3,2', 'C1,5', 'C4,4', 'C1,2', 'C2,2']
obstacles = ['C4,3', 'C3,5', 'C1,1', 'C2,3', 'C5,3', 'C2,1', 'C4,2', 'C5,1']

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

# Find the path
path = find_complete_path(start, goals, adjacency, obstacles)

# Output the result
print(f"<<<{json.dumps(path)}>>>")