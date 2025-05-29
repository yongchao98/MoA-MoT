from collections import deque
import heapq
import json

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def get_path(start, end, adjacency, obstacles):
    if start == end:
        return [start]
    
    # Priority queue for A* search
    pq = [(0, start, [start])]
    visited = set()
    
    while pq:
        _, current, path = heapq.heappop(pq)
        
        if current == end:
            return path
            
        if current in visited:
            continue
            
        visited.add(current)
        
        for next_pos in adjacency[current]:
            if next_pos not in obstacles and next_pos not in visited:
                new_path = path + [next_pos]
                priority = len(new_path) + manhattan_distance(next_pos, end)
                heapq.heappush(pq, (priority, next_pos, new_path))
    
    return None

def find_complete_path(start, goals, adjacency, obstacles):
    current = start
    unvisited_goals = set(goals)
    complete_path = [start]
    
    while unvisited_goals:
        # Find nearest unvisited goal
        min_dist = float('inf')
        nearest_goal = None
        
        for goal in unvisited_goals:
            dist = manhattan_distance(current, goal)
            if dist < min_dist:
                min_dist = dist
                nearest_goal = goal
        
        # Find path to nearest goal
        path = get_path(current, nearest_goal, adjacency, obstacles)
        if path is None:
            return None
        
        # Add path (excluding start position which is already included)
        complete_path.extend(path[1:])
        current = nearest_goal
        unvisited_goals.remove(nearest_goal)
    
    return complete_path

# Problem setup
start = "C1,5"
goals = ['C2,6', 'C5,6', 'C6,1', 'C2,1', 'C4,1', 'C3,5', 'C5,5']
obstacles = ['C3,1', 'C2,2', 'C4,5', 'C5,4', 'C5,2', 'C6,2', 'C1,6', 'C6,4', 'C4,6', 'C2,7']

# Adjacency dictionary
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C2,7"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"],
    "C2,6": ["C2,5", "C2,7", "C1,6", "C3,6"],
    "C2,7": ["C2,6", "C1,7", "C3,7"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"],
    "C3,6": ["C3,5", "C3,7", "C2,6", "C4,6"],
    "C3,7": ["C3,6", "C2,7", "C4,7"],
    "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5", "C5,5"],
    "C4,6": ["C4,5", "C4,7", "C3,6", "C5,6"],
    "C4,7": ["C4,6", "C3,7", "C5,7"],
    "C5,1": ["C5,2", "C4,1", "C6,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2", "C6,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3", "C6,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4", "C6,4"],
    "C5,5": ["C5,4", "C5,6", "C4,5", "C6,5"],
    "C5,6": ["C5,5", "C5,7", "C4,6", "C6,6"],
    "C5,7": ["C5,6", "C4,7", "C6,7"],
    "C6,1": ["C6,2", "C5,1"],
    "C6,2": ["C6,1", "C6,3", "C5,2"],
    "C6,3": ["C6,2", "C6,4", "C5,3"],
    "C6,4": ["C6,3", "C6,5", "C5,4"],
    "C6,5": ["C6,4", "C6,6", "C5,5"],
    "C6,6": ["C6,5", "C6,7", "C5,6"],
    "C6,7": ["C6,6", "C5,7"]
}

# Find solution
solution = find_complete_path(start, goals, adjacency, obstacles)

# Format and print solution
if solution:
    print(f"<<<{json.dumps(solution)}>>>")
else:
    print("No solution found")