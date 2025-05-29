import json
from heapq import heappush, heappop
from typing import List, Set, Dict, Tuple

def manhattan_distance(pos1: str, pos2: str) -> int:
    x1, y1 = map(int, pos1.split(','))
    x2, y2 = map(int, pos2.split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def get_coordinates(cell: str) -> Tuple[int, int]:
    return tuple(map(int, cell.replace('C', '').split(',')))

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
    unvisited_goals = set(goals)
    current_pos = start
    complete_path = [start]
    
    while unvisited_goals:
        # Find nearest unvisited goal
        min_dist = float('inf')
        nearest_goal = None
        
        for goal in unvisited_goals:
            path = a_star(current_pos, goal, adjacency, obstacles)
            if path and len(path) < min_dist:
                min_dist = len(path)
                nearest_goal = goal
        
        if nearest_goal is None:
            return []  # No valid path found
            
        # Add path to nearest goal (excluding start position to avoid duplication)
        path_to_goal = a_star(current_pos, nearest_goal, adjacency, obstacles)
        complete_path.extend(path_to_goal[1:])
        
        current_pos = nearest_goal
        unvisited_goals.remove(nearest_goal)
    
    return complete_path

# Main execution
start = "C5,4"
goals = ['C3,2', 'C6,3', 'C5,2', 'C3,6', 'C4,7', 'C3,7', 'C6,1']
obstacles = {'C6,7', 'C4,3', 'C6,2', 'C5,5', 'C3,1', 'C1,7', 'C2,2', 'C1,2', 'C6,4', 'C4,6'}

# Load adjacency data
adjacency = {
    # ... (your provided adjacency dict)
}

# Find the complete path
path = find_complete_path(start, goals, adjacency, obstacles)

# Format and print the result
result = json.dumps(path)
print(f"<<<{result}>>>")