import json
from heapq import heappush, heappop
from typing import List, Dict, Set, Tuple

def manhattan_distance(pos1: str, pos2: str) -> int:
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def format_cell(cell: str) -> str:
    if not cell.startswith('C'):
        return f"C{cell}"
    return cell

def a_star(start: str, goal: str, adjacency: Dict[str, List[str]], obstacles: Set[str]) -> List[str]:
    start = format_cell(start)
    goal = format_cell(goal)
    
    if start in obstacles or goal in obstacles:
        return []
    
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
    current = format_cell(start)
    goals = [format_cell(g) for g in goals]
    obstacles = {format_cell(o) for o in obstacles}
    
    unvisited_goals = set(goals)
    path = [current]
    
    while unvisited_goals:
        # Find nearest unvisited goal
        min_dist = float('inf')
        nearest_goal = None
        best_segment = None
        
        for goal in unvisited_goals:
            segment = a_star(current, goal, adjacency, obstacles)
            if segment and len(segment) < min_dist:
                min_dist = len(segment)
                nearest_goal = goal
                best_segment = segment
        
        if nearest_goal is None:
            return []  # No valid path found
        
        # Add the path segment (excluding the start position which is already in the path)
        path.extend(best_segment[1:])
        current = nearest_goal
        unvisited_goals.remove(nearest_goal)
    
    return path

# Main execution
start = "C2,6"
goals = ['C6,3', 'C1,4', 'C3,4', 'C2,2', 'C1,3', 'C4,5', 'C4,2']
obstacles = {'C3,3', 'C6,6', 'C6,2', 'C5,2', 'C1,5', 'C1,1', 'C2,1', 'C2,5', 'C3,2', 'C4,4'}

# Define adjacency dictionary
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
    "C5,1": ["C5,2", "C4,1", "C6,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2", "C6,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3", "C6,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4", "C6,4"],
    "C5,5": ["C5,4", "C5,6", "C4,5", "C6,5"],
    "C5,6": ["C5,5", "C4,6", "C6,6"],
    "C6,1": ["C6,2", "C5,1"],
    "C6,2": ["C6,1", "C6,3", "C5,2"],
    "C6,3": ["C6,2", "C6,4", "C5,3"],
    "C6,4": ["C6,3", "C6,5", "C5,4"],
    "C6,5": ["C6,4", "C6,6", "C5,5"],
    "C6,6": ["C6,5", "C5,6"]
}

# Find the complete path
final_path = find_complete_path(start, goals, adjacency, obstacles)

# Format the output
print(f"<<<{json.dumps(final_path)}>>>")