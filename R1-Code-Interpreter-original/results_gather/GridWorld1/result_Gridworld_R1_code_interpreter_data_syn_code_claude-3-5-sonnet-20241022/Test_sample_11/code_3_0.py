import json
from heapq import heappush, heappop
from typing import List, Dict, Set
from itertools import permutations

def manhattan_distance(pos1: str, pos2: str) -> int:
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def find_path(start: str, goal: str, adjacency: Dict, obstacles: Set[str], visited_goals: Set[str]) -> List[str]:
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
        
        # Sort adjacent positions by their distance to goal
        next_positions = sorted(
            [pos for pos in adjacency[current] if pos not in obstacles],
            key=lambda x: manhattan_distance(x, goal)
        )
        
        for next_pos in next_positions:
            if next_pos not in visited or next_pos in visited_goals:
                new_path = path + [next_pos]
                priority = len(new_path) + manhattan_distance(next_pos, goal)
                heappush(queue, (priority, next_pos, new_path))
    
    return []

def simplify_path(path: List[str], goals: Set[str]) -> List[str]:
    if len(path) <= 2:
        return path
        
    result = [path[0]]
    for i in range(1, len(path)-1):
        current = path[i]
        if current in goals or current == path[-1]:
            result.append(current)
        elif path[i-1] != path[i+1]:  # Keep points that are part of actual turns
            result.append(current)
    result.append(path[-1])
    return result

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

# Start with closest goals first
current_pos = start_pos
unvisited_goals = sorted(goals, key=lambda g: manhattan_distance(start_pos, g))
final_path = [start_pos]
visited_goals = set()

while unvisited_goals:
    goal = unvisited_goals[0]
    path = find_path(current_pos, goal, adjacency, obstacles, set(goals))
    
    if path:
        final_path.extend(path[1:])
        current_pos = goal
        visited_goals.add(goal)
        unvisited_goals.remove(goal)
    else:
        unvisited_goals.remove(goal)

# Simplify the path while maintaining goal visits
final_path = simplify_path(final_path, set(goals))

# Verify path validity
valid = True
for i in range(len(final_path)-1):
    if final_path[i+1] not in adjacency[final_path[i]] or final_path[i] in obstacles:
        valid = False
        break

if valid and all(goal in final_path for goal in goals):
    print(f"<<<{json.dumps(final_path)}>>>")
else:
    print("No valid path found")