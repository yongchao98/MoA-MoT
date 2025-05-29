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
    unvisited_goals = set(goals)
    path = [start]
    
    # Keep track of visited goals
    visited_goals = set()
    
    while unvisited_goals:
        # Find the nearest unvisited goal with minimum total path cost
        min_cost = float('inf')
        best_goal = None
        best_segment = None
        
        for goal in unvisited_goals:
            segment = find_path(current, goal, adjacency, obstacles)
            if segment and len(segment) < min_cost:
                min_cost = len(segment)
                best_goal = goal
                best_segment = segment
        
        if not best_segment:
            return []  # No valid path found
        
        # Add the path segment (excluding the start position which is already in the path)
        path.extend(best_segment[1:])
        current = best_goal
        unvisited_goals.remove(best_goal)
        visited_goals.add(best_goal)
        
        # If this is the last goal, don't add any more moves
        if not unvisited_goals:
            break
    
    return path

# Problem setup
initial_pos = "C5,1"
goals = ["C2,5", "C2,6", "C2,1", "C4,2", "C1,5", "C1,4", "C5,2"]
obstacles = {"C2,2", "C1,6", "C1,3", "C5,6", "C5,3", "C2,4", "C2,3", "C3,5", "C1,1", "C3,3"}

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

result = find_complete_path(initial_pos, goals, adjacency, obstacles)
print(f"<<<{json.dumps(result)}>>>")