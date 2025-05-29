import heapq
from typing import List, Dict, Set, Tuple
import json

def manhattan_distance(pos1: str, pos2: str) -> int:
    # Remove 'C' prefix and split coordinates
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def get_path_astar(start: str, goal: str, adjacency: Dict, obstacles: Set[str]) -> List[str]:
    if start == goal:
        return [start]
    
    frontier = [(0, start, [start])]
    visited = set()
    
    while frontier:
        _, current, path = heapq.heappop(frontier)
        
        if current == goal:
            return path
            
        if current in visited:
            continue
            
        visited.add(current)
        
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                new_path = path + [next_pos]
                priority = len(new_path) + manhattan_distance(next_pos, goal)
                heapq.heappush(frontier, (priority, next_pos, new_path))
    
    return []

def find_complete_path(start: str, goals: List[str], adjacency: Dict, obstacles: List[str]) -> List[str]:
    obstacles_set = set(obstacles)
    unvisited_goals = set(goals)
    current_pos = start
    complete_path = [start]
    
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
        path_segment = get_path_astar(current_pos, nearest_goal, adjacency, obstacles_set)
        
        if not path_segment:
            continue
            
        # Add path segment (excluding start position if not first segment)
        complete_path.extend(path_segment[1:])
        current_pos = nearest_goal
        unvisited_goals.remove(nearest_goal)
    
    return complete_path

# Problem inputs
start = "C6,6"
goals = ['C1,1', 'C3,4', 'C4,1', 'C4,4', 'C3,2', 'C1,4', 'C2,4']
obstacles = ['C2,5', 'C5,4', 'C3,1', 'C6,1', 'C4,6', 'C6,4', 'C3,3', 'C2,2', 'C1,6', 'C3,6']
adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"], "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"], "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C2,6"], "C2,1": ["C2,2", "C1,1", "C3,1"], "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"], "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"], "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"], "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"], "C2,6": ["C2,5", "C1,6", "C3,6"], "C3,1": ["C3,2", "C2,1", "C4,1"], "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"], "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"], "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"], "C3,6": ["C3,5", "C2,6", "C4,6"], "C4,1": ["C4,2", "C3,1", "C5,1"], "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"], "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"], "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"], "C4,5": ["C4,4", "C4,6", "C3,5", "C5,5"], "C4,6": ["C4,5", "C3,6", "C5,6"], "C5,1": ["C5,2", "C4,1", "C6,1"], "C5,2": ["C5,1", "C5,3", "C4,2", "C6,2"], "C5,3": ["C5,2", "C5,4", "C4,3", "C6,3"], "C5,4": ["C5,3", "C5,5", "C4,4", "C6,4"], "C5,5": ["C5,4", "C5,6", "C4,5", "C6,5"], "C5,6": ["C5,5", "C4,6", "C6,6"], "C6,1": ["C6,2", "C5,1"], "C6,2": ["C6,1", "C6,3", "C5,2"], "C6,3": ["C6,2", "C6,4", "C5,3"], "C6,4": ["C6,3", "C6,5", "C5,4"], "C6,5": ["C6,4", "C6,6", "C5,5"], "C6,6": ["C6,5", "C5,6"]
}

# Find solution
path = find_complete_path(start, goals, adjacency, obstacles)

# Output result
print(f"<<<{json.dumps(path)}>>>")