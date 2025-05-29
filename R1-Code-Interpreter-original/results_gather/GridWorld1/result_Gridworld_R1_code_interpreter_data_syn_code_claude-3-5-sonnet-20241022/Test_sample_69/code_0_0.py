import json
from heapq import heappush, heappop
from typing import List, Set, Dict, Tuple
import math

def manhattan_distance(pos1: str, pos2: str) -> int:
    x1, y1 = map(int, pos1.split(','))
    x2, y2 = map(int, pos2.split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def get_path(came_from: Dict[str, str], current: str) -> List[str]:
    path = [current]
    while current in came_from:
        current = came_from[current]
        path.append(current)
    return path[::-1]

def a_star(start: str, goal: str, adjacency: Dict[str, List[str]], obstacles: Set[str]) -> List[str]:
    if start == goal:
        return [start]
    
    frontier = [(0, start)]
    came_from = {}
    cost_so_far = {start: 0}
    
    while frontier:
        _, current = heappop(frontier)
        
        if current == goal:
            return get_path(came_from, goal)
            
        for next_pos in adjacency[current]:
            if next_pos in obstacles:
                continue
                
            new_cost = cost_so_far[current] + 1
            if next_pos not in cost_so_far or new_cost < cost_so_far[next_pos]:
                cost_so_far[next_pos] = new_cost
                priority = new_cost + manhattan_distance(next_pos, goal)
                heappush(frontier, (priority, next_pos))
                came_from[next_pos] = current
    
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
            
        # Add path to nearest goal to complete path
        path = a_star(current_pos, nearest_goal, adjacency, obstacles)
        complete_path.extend(path[1:])  # Skip first position as it's already in the path
        
        current_pos = nearest_goal
        unvisited_goals.remove(nearest_goal)
    
    return complete_path

# Main execution
start = "C3,2"
goals = ["C2,4", "C2,5", "C1,6", "C2,3", "C6,6", "C6,3", "C6,2"]
obstacles = {"C1,7", "C2,7", "C6,7", "C5,6", "C3,5", "C4,3", "C1,1", "C2,6", "C4,7", "C2,2"}

# Load adjacency data
adjacency = {
    # ... (using the provided adjacency dictionary)
}

# Find and print the path
path = find_complete_path(start, goals, adjacency, obstacles)
print(f"<<<{json.dumps(path)}>>>")