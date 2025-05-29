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
    
    # Priority queue for A* search
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

# Main problem data
start_pos = "C3,2"
goals = ['C1,4', 'C1,2', 'C3,5', 'C2,5', 'C2,1', 'C3,3', 'C2,2']
obstacles = {'C4,5', 'C4,6', 'C5,4', 'C2,3', 'C1,1', 'C5,6', 'C3,1', 'C5,3', 'C5,2', 'C4,3'}

# Load adjacency data
adjacency = {
    # ... (your provided adjacency dict)
}

# Find complete path visiting all goals
current_pos = start_pos
unvisited_goals = goals.copy()
final_path = [start_pos]
visited_goals = set()

while unvisited_goals:
    # Find nearest unvisited goal
    nearest_goal = min(unvisited_goals, key=lambda g: manhattan_distance(current_pos, g))
    
    # Find path to nearest goal
    path = find_path(current_pos, nearest_goal, adjacency, obstacles)
    
    if path:
        # Add path (excluding start position which is already included)
        final_path.extend(path[1:])
        current_pos = nearest_goal
        unvisited_goals.remove(nearest_goal)
        visited_goals.add(nearest_goal)
    else:
        # If no path found to this goal, try next closest goal
        unvisited_goals.remove(nearest_goal)

# Verify path validity
valid = True
for i in range(len(final_path)-1):
    if final_path[i+1] not in adjacency[final_path[i]]:
        valid = False
        break
    if final_path[i] in obstacles:
        valid = False
        break

if valid and len(visited_goals) == len(goals):
    print(f"<<<{json.dumps(final_path)}>>>")
else:
    print("No valid path found")