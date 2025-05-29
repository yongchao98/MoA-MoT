import json
from collections import deque
from itertools import permutations

def find_path(start, end, adjacency, obstacles):
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

# Define the problem
start = "C3,3"
goals = ["C1,1", "C4,1", "C3,2", "C1,3"]
obstacles = ["C3,5", "C2,4", "C4,3", "C2,3", "C4,5"]

# Define adjacency (using the provided adjacency dict)
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
    "C4,1": ["C4,2", "C3,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C4,5": ["C4,4", "C3,5"]
}

# Try all possible goal orderings
best_path = None
min_length = float('inf')

for goal_order in permutations(goals):
    current_path = [start]
    current_pos = start
    valid_sequence = True
    
    for goal in goal_order:
        path_segment = find_path(current_pos, goal, adjacency, obstacles)
        if not path_segment:
            valid_sequence = False
            break
        current_path.extend(path_segment[1:])  # Skip first position as it's the same as current_pos
        current_pos = goal
    
    if valid_sequence and len(current_path) < min_length:
        min_length = len(current_path)
        best_path = current_path

print(f"<<<{json.dumps(best_path)}>>>")