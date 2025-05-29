import json
from heapq import heappush, heappop

def manhattan_distance(pos1: str, pos2: str) -> int:
    # Remove the 'C' prefix and split coordinates
    x1, y1 = map(int, pos1.replace('C', '').split(','))
    x2, y2 = map(int, pos2.replace('C', '').split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def get_path_astar(start: str, goal: str, adjacency: dict, obstacles: set) -> list:
    if start == goal:
        return [start]
    
    # Priority queue: (f_score, current_pos, path)
    open_set = [(manhattan_distance(start, goal), start, [start])]
    closed_set = set()
    
    while open_set:
        _, current, path = heappop(open_set)
        
        if current == goal:
            return path
            
        if current in closed_set:
            continue
            
        closed_set.add(current)
        
        for neighbor in adjacency[current]:
            if neighbor not in obstacles and neighbor not in closed_set:
                new_path = path + [neighbor]
                g_score = len(new_path) - 1
                h_score = manhattan_distance(neighbor, goal)
                f_score = g_score + h_score
                heappush(open_set, (f_score, neighbor, new_path))
    
    return []

def find_nearest_goal(current: str, goals: set, adjacency: dict, obstacles: set) -> tuple:
    min_dist = float('inf')
    nearest_goal = None
    best_path = None
    
    for goal in goals:
        path = get_path_astar(current, goal, adjacency, obstacles)
        if path and len(path) < min_dist:
            min_dist = len(path)
            nearest_goal = goal
            best_path = path
            
    return nearest_goal, best_path

def solve_gridworld(start: str, goals: list, adjacency: dict, obstacles: set) -> list:
    current = start
    unvisited_goals = set(goals)
    complete_path = [start]
    
    while unvisited_goals:
        next_goal, path = find_nearest_goal(current, unvisited_goals, adjacency, obstacles)
        
        if not path:
            return []  # No valid solution found
            
        # Add path to next goal (excluding the start position)
        complete_path.extend(path[1:])
        current = next_goal
        unvisited_goals.remove(next_goal)
    
    return complete_path

# Problem setup
start = "C2,6"
goals = ['C5,1', 'C2,3', 'C3,5', 'C5,6', 'C4,4', 'C5,2', 'C1,6']
obstacles = {'C5,3', 'C4,6', 'C4,1', 'C3,1', 'C1,3', 'C3,3', 'C1,4', 'C1,5', 'C2,4', 'C2,1'}

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

# Find solution
path = solve_gridworld(start, goals, adjacency, obstacles)

# Output result
print(f"<<<{json.dumps(path)}>>>")