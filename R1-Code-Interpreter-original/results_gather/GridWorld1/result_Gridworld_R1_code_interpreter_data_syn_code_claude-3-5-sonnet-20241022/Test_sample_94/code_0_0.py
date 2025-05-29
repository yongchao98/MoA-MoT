import json
from collections import deque
import itertools

# Grid setup
obstacles = ['C4,2', 'C2,2', 'C3,2', 'C5,2', 'C4,3', 'C3,3', 'C5,3', 'C2,3']
goals = ['C1,4', 'C4,4', 'C2,1', 'C3,1', 'C3,4', 'C3,5']
start = 'C1,2'

# Adjacency list (given in the problem)
adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C2,5"], "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"], "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"], "C2,5": ["C2,4", "C1,5", "C3,5"],
    "C3,1": ["C3,2", "C2,1", "C4,1"], "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C2,5", "C4,5"], "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"], "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"], "C4,5": ["C4,4", "C3,5", "C5,5"],
    "C5,1": ["C5,2", "C4,1"], "C5,2": ["C5,1", "C5,3", "C4,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3"], "C5,4": ["C5,3", "C5,5", "C4,4"],
    "C5,5": ["C5,4", "C4,5"]
}

def find_path(start, end, obstacles):
    if start == end:
        return [start]
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                if next_pos == end:
                    return path + [next_pos]
                visited.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def find_complete_path():
    # Find all possible paths between goals
    paths_between_goals = {}
    all_points = [start] + goals
    
    for point1, point2 in itertools.combinations(all_points, 2):
        path = find_path(point1, point2, obstacles)
        if path:
            paths_between_goals[(point1, point2)] = path
            paths_between_goals[(point2, point1)] = path[::-1]
    
    # Try different goal orderings
    best_path = None
    min_length = float('inf')
    
    for goal_order in itertools.permutations(goals):
        current_path = [start]
        current_pos = start
        valid_sequence = True
        
        for goal in goal_order:
            if (current_pos, goal) not in paths_between_goals:
                valid_sequence = False
                break
            path_segment = paths_between_goals[(current_pos, goal)][1:]  # Exclude start position
            current_path.extend(path_segment)
            current_pos = goal
            
        if valid_sequence and len(current_path) < min_length:
            min_length = len(current_path)
            best_path = current_path
    
    return best_path

# Find and print the solution
solution = find_complete_path()
print(f"<<<{json.dumps(solution)}>>>")