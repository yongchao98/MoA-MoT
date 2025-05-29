import json
from heapq import heappush, heappop
from typing import List, Dict, Set, Tuple

def manhattan_distance(pos1: str, pos2: str) -> int:
    x1, y1 = map(int, pos1.split(','))
    x2, y2 = map(int, pos2.split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def get_path_astar(start: str, goal: str, adjacency: Dict, obstacles: Set[str]) -> List[str]:
    if start == goal:
        return [start]
    
    open_set = [(0, start, [start])]
    closed_set = set()
    
    while open_set:
        f_score, current, path = heappop(open_set)
        
        if current == goal:
            return path
            
        if current in closed_set:
            continue
            
        closed_set.add(current)
        
        for neighbor in adjacency[current]:
            if neighbor not in obstacles and neighbor not in closed_set:
                new_path = path + [neighbor]
                h_score = manhattan_distance(neighbor, goal)
                g_score = len(new_path) - 1
                f_score = g_score + h_score
                heappush(open_set, (f_score, neighbor, new_path))
    
    return []

def find_complete_path(start: str, goals: List[str], adjacency: Dict, obstacles: Set[str]) -> List[str]:
    current = start
    unvisited_goals = set(goals)
    complete_path = [start]
    
    while unvisited_goals:
        # Find nearest unvisited goal
        min_dist = float('inf')
        next_goal = None
        next_path = None
        
        for goal in unvisited_goals:
            path = get_path_astar(current, goal, adjacency, obstacles)
            if path and len(path) < min_dist:
                min_dist = len(path)
                next_goal = goal
                next_path = path
        
        if next_goal is None:
            return []  # No valid path found
            
        # Add path to next goal (excluding the start position which is already in the path)
        complete_path.extend(next_path[1:])
        current = next_goal
        unvisited_goals.remove(next_goal)
    
    return complete_path

# Main execution
start = "C2,6"
goals = ['C5,1', 'C2,3', 'C3,5', 'C5,6', 'C4,4', 'C5,2', 'C1,6']
obstacles = {'C5,3', 'C4,6', 'C4,1', 'C3,1', 'C1,3', 'C3,3', 'C1,4', 'C1,5', 'C2,4', 'C2,1'}

# Load adjacency data
adjacency = {
    # ... (using the provided adjacency dictionary)
}

# Find path
path = find_complete_path(start, goals, adjacency, obstacles)

# Output result
print(f"<<<{json.dumps(path)}>>>")