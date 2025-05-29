import json
from heapq import heappush, heappop
from typing import List, Dict, Set, Tuple

def manhattan_distance(pos1: str, pos2: str) -> int:
    x1, y1 = map(int, pos1.split(','))
    x2, y2 = map(int, pos2.split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def get_coordinates(pos: str) -> Tuple[int, int]:
    return tuple(map(int, pos[1:].split(',')))

def a_star(start: str, goal: str, adjacency: Dict[str, List[str]], obstacles: Set[str]) -> List[str]:
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

def find_complete_path(start: str, goals: List[str], adjacency: Dict[str, List[str]], 
                      obstacles: List[str]) -> List[str]:
    obstacles_set = set(obstacles)
    unvisited_goals = set(goals)
    current_pos = start
    complete_path = [start]
    
    while unvisited_goals:
        # Find nearest unvisited goal
        min_dist = float('inf')
        next_goal = None
        
        for goal in unvisited_goals:
            path = a_star(current_pos, goal, adjacency, obstacles_set)
            if path and len(path) < min_dist:
                min_dist = len(path)
                next_goal = goal
        
        if next_goal is None:
            return []  # No valid path found
            
        # Add path to next goal (excluding start position to avoid duplicates)
        path_to_goal = a_star(current_pos, next_goal, adjacency, obstacles_set)
        complete_path.extend(path_to_goal[1:])
        
        current_pos = next_goal
        unvisited_goals.remove(next_goal)
    
    return complete_path

# Main execution
start = "C3,1"
goals = ['C1,1', 'C1,7', 'C4,6', 'C2,3', 'C1,5', 'C4,7', 'C3,7']
obstacles = ['C2,4', 'C3,6', 'C4,4', 'C2,7', 'C2,2', 'C5,6', 'C4,3', 'C4,2', 'C5,3', 'C5,2']

# Load adjacency data
adjacency = {
    # ... (your provided adjacency dict)
}

# Find and print the path
path = find_complete_path(start, goals, adjacency, obstacles)
print(f"<<<{json.dumps(path)}>>>")