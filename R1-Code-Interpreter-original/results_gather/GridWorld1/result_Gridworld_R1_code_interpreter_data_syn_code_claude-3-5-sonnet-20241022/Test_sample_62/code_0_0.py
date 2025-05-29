import heapq
from typing import List, Dict, Set, Tuple
import json

def manhattan_distance(pos1: str, pos2: str) -> int:
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def a_star(start: str, goal: str, adjacency: Dict, obstacles: Set[str]) -> List[str]:
    if start in obstacles or goal in obstacles:
        return []
    
    frontier = [(0, start, [start])]
    visited = set()
    
    while frontier:
        f_score, current, path = heapq.heappop(frontier)
        
        if current == goal:
            return path
            
        if current in visited:
            continue
            
        visited.add(current)
        
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                new_path = path + [next_pos]
                h_score = manhattan_distance(next_pos, goal)
                g_score = len(new_path)
                f_score = g_score + h_score
                heapq.heappush(frontier, (f_score, next_pos, new_path))
    
    return []

def find_complete_path(start: str, goals: List[str], adjacency: Dict, obstacles: Set[str]) -> List[str]:
    current_pos = start
    unvisited_goals = goals.copy()
    final_path = [start]
    
    while unvisited_goals:
        # Find nearest unvisited goal
        min_dist = float('inf')
        nearest_goal = None
        best_path = None
        
        for goal in unvisited_goals:
            path = a_star(current_pos, goal, adjacency, obstacles)
            if path and len(path) < min_dist:
                min_dist = len(path)
                nearest_goal = goal
                best_path = path
        
        if nearest_goal is None:
            return []  # No valid path found
        
        # Add path to nearest goal (excluding start position if not first path)
        final_path.extend(best_path[1:])
        current_pos = nearest_goal
        unvisited_goals.remove(nearest_goal)
    
    return final_path

# Problem setup
start = "C4,1"
goals = ['C6,4', 'C3,2', 'C5,5', 'C3,7', 'C6,7', 'C2,7', 'C4,7']
obstacles = {'C1,1', 'C6,1', 'C1,2', 'C4,6', 'C3,3', 'C6,6', 'C2,1', 'C1,5', 'C2,3', 'C3,5'}

# Load adjacency data
adjacency = {
    # ... (using the provided adjacency dictionary)
}

# Find solution
solution = find_complete_path(start, goals, adjacency, obstacles)

# Format output
print(f"<<<{json.dumps(solution)}>>>")