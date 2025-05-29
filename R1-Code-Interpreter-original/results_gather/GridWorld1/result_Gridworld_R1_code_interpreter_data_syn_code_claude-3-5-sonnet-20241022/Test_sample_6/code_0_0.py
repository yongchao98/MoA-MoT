import json
from heapq import heappush, heappop
from typing import List, Dict, Set, Tuple
import math

def manhattan_distance(pos1: str, pos2: str) -> int:
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def find_path(start: str, goal: str, adjacency: Dict, obstacles: Set[str]) -> List[str]:
    if start == goal:
        return [start]
    
    # A* algorithm
    frontier = [(0, start, [start])]
    visited = set()
    
    while frontier:
        _, current, path = heappop(frontier)
        
        if current == goal:
            return path
            
        if current in visited:
            continue
            
        visited.add(current)
        
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                new_path = path + [next_pos]
                priority = len(new_path) + manhattan_distance(next_pos, goal)
                heappush(frontier, (priority, next_pos, new_path))
    
    return []

def find_complete_path(start: str, goals: List[str], adjacency: Dict, obstacles: Set[str]) -> List[str]:
    current = start
    unvisited_goals = set(goals)
    path = [start]
    
    while unvisited_goals:
        # Find nearest unvisited goal
        min_dist = float('inf')
        nearest_goal = None
        
        for goal in unvisited_goals:
            dist = manhattan_distance(current, goal)
            if dist < min_dist:
                min_dist = dist
                nearest_goal = goal
        
        if nearest_goal:
            # Find path to nearest goal
            segment = find_path(current, nearest_goal, adjacency, obstacles)
            if segment:
                path.extend(segment[1:])  # Skip first position as it's already in path
                current = nearest_goal
                unvisited_goals.remove(nearest_goal)
            else:
                # If no path found, try next nearest goal
                unvisited_goals.remove(nearest_goal)
                continue
    
    return path

# Main execution
start = "C3,3"
goals = ['C1,1', 'C2,1', 'C5,5', 'C3,5', 'C6,4', 'C6,3', 'C5,1']
obstacles = {'C5,3', 'C1,2', 'C4,2', 'C2,7', 'C3,7', 'C2,2', 'C1,4', 'C6,2', 'C5,6', 'C4,5'}

# Load adjacency data
adjacency = {
    # ... (using the provided adjacency dictionary)
}

# Find complete path
final_path = find_complete_path(start, goals, adjacency, obstacles)

# Output result
print(f"<<<{json.dumps(final_path)}>>>")