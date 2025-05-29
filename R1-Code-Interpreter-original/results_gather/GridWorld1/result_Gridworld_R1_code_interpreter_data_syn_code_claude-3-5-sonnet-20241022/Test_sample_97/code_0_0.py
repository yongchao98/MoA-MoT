import json
from heapq import heappush, heappop
from typing import List, Dict, Set, Tuple

def manhattan_distance(pos1: str, pos2: str) -> int:
    x1, y1 = map(int, pos1.split(','))
    x2, y2 = map(int, pos2.split(','))
    return abs(x1 - y1) + abs(x2 - y2)

def get_coordinates(cell: str) -> Tuple[int, int]:
    # Remove 'C' prefix and split coordinates
    return tuple(map(int, cell[1:].split(',')))

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
                      obstacles: Set[str]) -> List[str]:
    current = start
    unvisited_goals = set(goals)
    path = [start]
    
    while unvisited_goals:
        # Find nearest unvisited goal
        min_dist = float('inf')
        next_goal = None
        next_path = None
        
        for goal in unvisited_goals:
            temp_path = a_star(current, goal, adjacency, obstacles)
            if temp_path and len(temp_path) < min_dist:
                min_dist = len(temp_path)
                next_goal = goal
                next_path = temp_path
        
        if next_goal is None:
            return []  # No valid path found
            
        # Add path to next goal (excluding start position which is already in path)
        path.extend(next_path[1:])
        current = next_goal
        unvisited_goals.remove(next_goal)
    
    return path

# Main execution
start = "C2,1"
goals = ['C2,7', 'C4,3', 'C3,3', 'C5,3', 'C4,4', 'C2,5', 'C6,5']
obstacles = {'C6,7', 'C5,4', 'C2,2', 'C3,6', 'C3,4', 'C6,1', 'C2,4', 'C4,2', 'C2,6', 'C5,7'}

# Load adjacency data
adjacency = {
    # ... (your provided adjacency dict)
}

# Find complete path
final_path = find_complete_path(start, goals, adjacency, obstacles)

# Format output
print(f"<<<{json.dumps(final_path)}>>>")