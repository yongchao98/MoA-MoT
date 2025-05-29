import json
from heapq import heappush, heappop
from typing import List, Dict, Set, Tuple

def manhattan_distance(pos1: str, pos2: str) -> int:
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def find_path(start: str, goal: str, adjacency: Dict, obstacles: Set[str]) -> List[str]:
    if start == goal:
        return [start]
    
    queue = [(0, 0, start, [start])]
    visited = set()
    
    while queue:
        _, cost, current, path = heappop(queue)
        
        if current == goal:
            return path
            
        if current in visited:
            continue
            
        visited.add(current)
        
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                new_path = path + [next_pos]
                new_cost = cost + 1
                priority = new_cost + manhattan_distance(next_pos, goal)
                heappush(queue, (priority, new_cost, next_pos, new_path))
    
    return []

def find_complete_path(start: str, goals: List[str], adjacency: Dict, obstacles: Set[str]) -> List[str]:
    current = start
    remaining_goals = set(goals)
    path = [start]
    
    while remaining_goals:
        nearest = min(remaining_goals, key=lambda g: manhattan_distance(current, g))
        segment = find_path(current, nearest, adjacency, obstacles)[1:]
        if not segment:
            return []
            
        path.extend(segment)
        current = nearest
        remaining_goals.remove(nearest)
    
    return path

# Problem setup
initial_pos = "C5,1"
goals = ["C2,5", "C2,6", "C2,1", "C4,2", "C1,5", "C1,4", "C5,2"]
obstacles = {"C2,2", "C1,6", "C1,3", "C5,6", "C5,3", "C2,4", "C2,3", "C3,5", "C1,1", "C3,3"}

# Complete adjacency dictionary
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"],
    "C2,6": ["C2,5", "C1,6", "C3,6"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"],
    "C3,6": ["C3,5", "C2,6", "C4,6"],
    "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5", "C5,5"],
    "C4,6": ["C4,5", "C3,6", "C5,6"],
    "C5,1": ["C5,2", "C4,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4"],
    "C5,5": ["C5,4", "C5,6", "C4,5"],
    "C5,6": ["C5,5", "C4,6"]
}

# Find and output the path
result = find_complete_path(initial_pos, goals, adjacency, obstacles)
print(f"<<<{json.dumps(result)}>>>")