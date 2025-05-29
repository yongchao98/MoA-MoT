import json
from collections import deque

def find_path(start, target, adjacency, obstacles, visited=None):
    if visited is None:
        visited = set()
    
    queue = deque([(start, [start])])
    seen = {start}
    
    while queue:
        current, path = queue.popleft()
        
        if current == target:
            return path
            
        for next_pos in adjacency[current]:
            if (next_pos not in seen and 
                next_pos not in obstacles and 
                (visited is None or next_pos not in visited)):
                seen.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    
    return None

def get_distance(pos1, pos2):
    x1, y1 = map(int, pos1.split(',')[0][1:]), map(int, pos1.split(',')[1])
    x2, y2 = map(int, pos2.split(',')[0][1:]), map(int, pos2.split(',')[1])
    return abs(list(x1)[0] - list(x2)[0]) + abs(list(y1)[0] - list(y2)[0])

def find_solution_path(start, goals, adjacency, obstacles):
    current = start
    path = [start]
    remaining_goals = set(goals)
    
    while remaining_goals:
        # Find the nearest accessible goal
        best_goal = None
        best_path = None
        min_total_cost = float('inf')
        
        for goal in remaining_goals:
            # Try to find a path to this goal
            temp_path = find_path(current, goal, adjacency, obstacles)
            if temp_path:
                # Calculate cost as path length
                cost = len(temp_path)
                if cost < min_total_cost:
                    min_total_cost = cost
                    best_path = temp_path
                    best_goal = goal
        
        if best_path is None:
            # Try alternative paths if direct path is not found
            for intermediate in adjacency[current]:
                if intermediate not in obstacles:
                    for goal in remaining_goals:
                        temp_path = find_path(intermediate, goal, adjacency, obstacles)
                        if temp_path:
                            full_path = [current] + temp_path
                            cost = len(full_path)
                            if cost < min_total_cost:
                                min_total_cost = cost
                                best_path = full_path
                                best_goal = goal
        
        if best_path is None:
            return None
        
        # Add the path (excluding the start position which is already in the path)
        path.extend(best_path[1:])
        current = best_goal
        remaining_goals.remove(best_goal)
    
    return path

# Problem setup
start = "C4,4"
goals = ['C4,2', 'C5,1', 'C5,4', 'C3,2', 'C3,3', 'C5,3']
obstacles = ['C1,2', 'C2,5', 'C5,2', 'C1,5', 'C3,5', 'C1,3', 'C1,4', 'C4,3']

adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C2,5"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C1,5", "C3,5"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C2,5", "C4,5"],
    "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C3,5", "C5,5"],
    "C5,1": ["C5,2", "C4,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4"],
    "C5,5": ["C5,4", "C4,5"]
}

solution = find_solution_path(start, goals, adjacency, obstacles)
if solution:
    print(f"<<<{json.dumps(solution)}>>>")
else:
    print("No valid path found")