import json
from collections import deque
import heapq

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def get_path_astar(start, goal, adjacency, obstacles):
    if start == goal:
        return [start]
    
    def h(pos):
        return manhattan_distance(pos, goal)
    
    open_set = [(0 + h(start), 0, start, [start])]
    closed_set = set()
    
    while open_set:
        f, g, current, path = heapq.heappop(open_set)
        
        if current == goal:
            return path
            
        if current in closed_set:
            continue
            
        closed_set.add(current)
        
        for neighbor in adjacency[current]:
            if neighbor not in obstacles and neighbor not in closed_set:
                new_g = g + 1
                new_path = path + [neighbor]
                heapq.heappush(open_set, (new_g + h(neighbor), new_g, neighbor, new_path))
    
    return None

def find_complete_path(start, goals, adjacency, obstacles):
    current_pos = start
    unvisited_goals = set(goals)
    final_path = [start]
    
    while unvisited_goals:
        # Find nearest unvisited goal
        min_dist = float('inf')
        nearest_goal = None
        
        for goal in unvisited_goals:
            dist = manhattan_distance(current_pos, goal)
            if dist < min_dist:
                min_dist = dist
                nearest_goal = goal
        
        # Find path to nearest goal
        path = get_path_astar(current_pos, nearest_goal, adjacency, obstacles)
        
        if path is None:
            continue
            
        # Add path (excluding start position) to final path
        final_path.extend(path[1:])
        current_pos = nearest_goal
        unvisited_goals.remove(nearest_goal)
    
    return final_path

# Define the problem parameters
start = "C3,4"
goals = ['C5,7', 'C1,5', 'C4,2', 'C2,1', 'C5,3', 'C3,7', 'C4,6']
obstacles = ['C6,1', 'C4,1', 'C5,2', 'C1,1', 'C1,7', 'C2,6', 'C5,1', 'C6,6', 'C1,6', 'C3,2']

# Load adjacency data
adjacency = {
    # ... (your provided adjacency dict)
}

# Find the complete path
path = find_complete_path(start, goals, adjacency, obstacles)

# Format the output
print(f"<<<{json.dumps(path)}>>>")