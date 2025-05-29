import json
from heapq import heappush, heappop
from typing import List, Set, Dict, Tuple

def manhattan_distance(pos1: str, pos2: str) -> int:
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def a_star(start: str, goal: str, adjacency: Dict, obstacles: Set[str]) -> List[str]:
    if start == goal:
        return [start]
    
    frontier = [(0, start, [start])]
    visited = set()
    
    while frontier:
        _, current, path = heappop(frontier)
        
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
    
    return []

def find_complete_path(start: str, goals: List[str], adjacency: Dict, obstacles: Set[str]) -> List[str]:
    unvisited_goals = set(goals)
    current_pos = start
    final_path = [start]
    
    while unvisited_goals:
        # Find nearest unvisited goal
        min_dist = float('inf')
        nearest_goal = None
        
        for goal in unvisited_goals:
            dist = manhattan_distance(current_pos, goal)
            if dist < min_dist:
                min_dist = dist
                nearest_goal = goal
        
        # Find path to nearest goal
        path = a_star(current_pos, nearest_goal, adjacency, obstacles)
        if not path:
            continue
            
        # Add path (excluding start position) to final path
        final_path.extend(path[1:])
        current_pos = nearest_goal
        unvisited_goals.remove(nearest_goal)
    
    return final_path

# Main execution
start = "C3,5"
goals = ['C4,6', 'C3,6', 'C5,5', 'C6,1', 'C5,1', 'C2,1', 'C1,2']
obstacles = {'C3,3', 'C1,3', 'C2,5', 'C6,5', 'C5,3', 'C4,7', 'C2,4', 'C4,2', 'C1,4', 'C5,7'}

# Load adjacency data
adjacency = {
    # ... (using the provided adjacency dictionary)
}

# Find complete path
path = find_complete_path(start, goals, adjacency, obstacles)

# Output result
print(f"<<<{json.dumps(path)}>>>")