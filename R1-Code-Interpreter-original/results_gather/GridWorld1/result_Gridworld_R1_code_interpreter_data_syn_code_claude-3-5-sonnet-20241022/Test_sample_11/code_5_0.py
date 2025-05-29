import json
from heapq import heappush, heappop
from typing import List, Dict, Set

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
        
        # Sort adjacent positions by distance to goal for better efficiency
        next_positions = sorted(
            [pos for pos in adjacency[current] if pos not in obstacles and pos not in visited],
            key=lambda x: manhattan_distance(x, goal)
        )
        
        for next_pos in next_positions:
            new_path = path + [next_pos]
            priority = len(new_path) + manhattan_distance(next_pos, goal)
            heappush(queue, (priority, next_pos, new_path))
    
    return []

start_pos = "C3,2"
goals = ['C3,3', 'C2,2', 'C2,1', 'C1,2', 'C1,4', 'C2,5', 'C3,5']  # Reordered for efficiency
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

# Find path visiting all goals in sequence
current_pos = start_pos
final_path = [start_pos]
visited_goals = set()

for goal in goals:
    path = find_path(current_pos, goal, adjacency, obstacles)
    if path:
        final_path.extend(path[1:])
        current_pos = goal
        visited_goals.add(goal)

# Verify path validity
valid = True
for i in range(len(final_path)-1):
    if final_path[i+1] not in adjacency[final_path[i]] or final_path[i] in obstacles:
        valid = False
        break

if valid and len(visited_goals) == len(goals):
    print(f"<<<{json.dumps(final_path)}>>>")
else:
    print("No valid path found")