import json
from heapq import heappush, heappop
from typing import List, Set, Dict, Tuple

def manhattan_distance(pos1: str, pos2: str) -> int:
    x1, y1 = map(int, pos1.split(','))
    x2, y2 = map(int, pos2.split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def get_coordinates(pos: str) -> Tuple[int, int]:
    return tuple(map(int, pos[1:].split(',')))

def a_star(start: str, goal: str, adjacency: Dict, obstacles: Set[str]) -> List[str]:
    if start == goal:
        return [start]
    
    frontier = [(0, start, [start])]
    visited = {start}
    
    while frontier:
        _, current, path = heappop(frontier)
        
        if current == goal:
            return path
            
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                visited.add(next_pos)
                new_path = path + [next_pos]
                priority = len(new_path) + manhattan_distance(next_pos, goal)
                heappush(frontier, (priority, next_pos, new_path))
    
    return []

def find_complete_path(start: str, goals: List[str], adjacency: Dict, obstacles: List[str]) -> List[str]:
    obstacles_set = set(obstacles)
    unvisited_goals = set(goals)
    current_pos = start
    complete_path = [start]
    
    while unvisited_goals:
        # Find nearest unvisited goal
        nearest_goal = min(unvisited_goals, key=lambda g: manhattan_distance(current_pos, g))
        
        # Find path to nearest goal
        path_segment = a_star(current_pos, nearest_goal, adjacency, obstacles_set)
        
        if not path_segment:
            # If no path found, try another goal
            continue
            
        # Add path segment (excluding start position to avoid duplication)
        complete_path.extend(path_segment[1:])
        current_pos = nearest_goal
        unvisited_goals.remove(nearest_goal)
    
    return complete_path

# Problem inputs
start = "C6,6"
goals = ['C4,6', 'C3,4', 'C1,6', 'C6,3', 'C2,4', 'C3,5', 'C2,3']
obstacles = ['C6,1', 'C3,2', 'C3,6', 'C4,4', 'C1,2', 'C6,5', 'C5,5', 'C1,3', 'C4,1', 'C5,1']

# Load adjacency data
adjacency = {
    # ... (your provided adjacency dict)
}

# Find and print solution
path = find_complete_path(start, goals, adjacency, obstacles)
print(f"<<<{json.dumps(path)}>>>")