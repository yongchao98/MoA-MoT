import json
from collections import deque
import heapq

# Helper function to get Manhattan distance between two cells
def manhattan_distance(cell1, cell2):
    x1, y1 = map(int, cell1.replace('C', '').split(','))
    x2, y2 = map(int, cell2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

# A* search between two points
def find_path(start, end, adjacency, obstacles):
    if start == end:
        return [start]
    
    heap = [(manhattan_distance(start, end), 0, [start])]
    visited = set()
    
    while heap:
        _, cost, path = heapq.heappop(heap)
        current = path[-1]
        
        if current == end:
            return path
            
        if current in visited:
            continue
            
        visited.add(current)
        
        for next_cell in adjacency[current]:
            if next_cell not in obstacles and next_cell not in visited:
                new_path = path + [next_cell]
                new_cost = cost + 1
                priority = new_cost + manhattan_distance(next_cell, end)
                heapq.heappush(heap, (priority, new_cost, new_path))
    
    return None

# Find nearest unvisited goal
def find_nearest_goal(current, unvisited_goals):
    return min(unvisited_goals, key=lambda g: manhattan_distance(current, g))

# Main solution
def solve_gridworld(start, goals, obstacles, adjacency):
    current = start
    unvisited_goals = set(goals)
    final_path = [start]
    
    while unvisited_goals:
        # Find nearest unvisited goal
        next_goal = find_nearest_goal(current, unvisited_goals)
        
        # Find path to this goal
        path = find_path(current, next_goal, adjacency, set(obstacles))
        
        if path is None:
            return None
        
        # Add path (excluding start position) to final path
        final_path.extend(path[1:])
        
        # Update current position and mark goal as visited
        current = next_goal
        unvisited_goals.remove(next_goal)
    
    return final_path

# Input data
start = "C3,5"
goals = ['C1,1', 'C4,3', 'C3,4', 'C2,1', 'C5,1', 'C1,3', 'C4,4']
obstacles = ['C3,6', 'C4,6', 'C3,1', 'C5,5', 'C5,3', 'C2,4', 'C2,3', 'C1,5', 'C5,4', 'C1,4']
adjacency = {
    # ... (using the provided adjacency dictionary)
}

# Solve and output
path = solve_gridworld(start, goals, obstacles, adjacency)
print(f"<<<{json.dumps(path)}>>>")