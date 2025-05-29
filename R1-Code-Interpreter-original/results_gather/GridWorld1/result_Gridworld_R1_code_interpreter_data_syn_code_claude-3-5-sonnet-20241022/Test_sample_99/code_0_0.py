import json
from collections import deque
from itertools import permutations

def find_path(start, end, adjacency, obstacles):
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        if current == end:
            return path
            
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                visited.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def find_complete_path(start, goals, adjacency, obstacles):
    # Try all possible goal orderings
    shortest_path = None
    min_length = float('inf')
    
    for goal_order in permutations(goals):
        current_path = [start]
        current_pos = start
        valid_sequence = True
        
        # Try to connect all goals in the current order
        for goal in goal_order:
            if current_pos != goal:
                segment = find_path(current_pos, goal, adjacency, obstacles)
                if segment is None:
                    valid_sequence = False
                    break
                current_path.extend(segment[1:])  # Exclude the first position as it's already included
                current_pos = goal
        
        if valid_sequence and len(current_path) < min_length:
            shortest_path = current_path
            min_length = len(current_path)
    
    return shortest_path

# Define the problem parameters
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

start = "C2,4"
goals = ["C2,2", "C4,3", "C3,3", "C4,1"]
obstacles = ["C2,1", "C2,3", "C4,4", "C3,1", "C3,5"]

# Find the path
path = find_complete_path(start, goals, adjacency, obstacles)

# Output the result
print(f"<<<{json.dumps(path)}>>>")