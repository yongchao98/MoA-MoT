import json
from heapq import heappush, heappop
from typing import List, Dict, Set
from itertools import permutations

def manhattan_distance(pos1: str, pos2: str) -> int:
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def find_path(start: str, goal: str, adjacency: Dict, obstacles: Set[str]) -> List[str]:
    if start == goal:
        return [start]
    
    queue = [(0, start, [start])]
    visited = set()
    
    while queue:
        _, current, path = heappop(queue)
        
        if current == goal:
            return path
            
        if current in visited:
            continue
            
        visited.add(current)
        
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                new_path = path + [next_pos]
                priority = len(new_path) + manhattan_distance(next_pos, goal)
                heappush(queue, (priority, next_pos, new_path))
    
    return []

def find_total_path_length(sequence, start, adjacency, obstacles):
    total_length = 0
    current = start
    path = [start]
    
    for goal in sequence:
        segment = find_path(current, goal, adjacency, obstacles)
        if not segment:
            return float('inf'), []
        total_length += len(segment) - 1
        path.extend(segment[1:])
        current = goal
    
    return total_length, path

start_pos = "C3,2"
goals = ['C1,4', 'C1,2', 'C3,5', 'C2,5', 'C2,1', 'C3,3', 'C2,2']
obstacles = {'C4,5', 'C4,6', 'C5,4', 'C2,3', 'C1,1', 'C5,6', 'C3,1', 'C5,3', 'C5,2', 'C4,3'}

adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"], "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"], "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C2,6"],
    "C2,1": ["C2,2", "C1,1", "C3,1"], "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"], "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"], "C2,6": ["C2,5", "C1,6", "C3,6"],
    "C3,1": ["C3,2", "C2,1", "C4,1"], "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"], "C3,6": ["C3,5", "C2,6", "C4,6"],
    "C4,1": ["C4,2", "C3,1", "C5,1"], "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"], "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5", "C5,5"], "C4,6": ["C4,5", "C3,6", "C5,6"],
    "C5,1": ["C5,2", "C4,1"], "C5,2": ["C5,1", "C5,3", "C4,2"], "C5,3": ["C5,2", "C5,4", "C4,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4"], "C5,5": ["C5,4", "C5,6", "C4,5"], "C5,6": ["C5,5", "C4,6"]
}

# Try different goal orderings to find the shortest path
best_length = float('inf')
best_path = None

# Try a subset of permutations to keep computation time reasonable
for goal_sequence in permutations(goals[:4]):
    remaining_goals = goals[4:]
    length, path = find_total_path_length(list(goal_sequence) + remaining_goals, start_pos, adjacency, obstacles)
    if length < best_length:
        best_length = length
        best_path = path

if best_path:
    print(f"<<<{json.dumps(best_path)}>>>")
else:
    print("No valid path found")