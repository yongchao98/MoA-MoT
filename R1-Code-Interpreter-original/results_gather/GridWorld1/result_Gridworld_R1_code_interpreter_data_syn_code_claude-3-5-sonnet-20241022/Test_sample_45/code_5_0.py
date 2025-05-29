import json
from collections import deque

# Initialize the problem
initial_pos = "C4,4"
goals = ['C4,2', 'C5,1', 'C5,4', 'C3,2', 'C3,3', 'C5,3']
obstacles = ['C1,2', 'C2,5', 'C5,2', 'C1,5', 'C3,5', 'C1,3', 'C1,4', 'C4,3']

# Adjacency dictionary
adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C2,5"], "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C1,5", "C3,5"], "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C2,5", "C4,5"], "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"],
    "C4,5": ["C4,4", "C3,5", "C5,5"], "C5,1": ["C5,2", "C4,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2"], "C5,3": ["C5,2", "C5,4", "C4,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4"], "C5,5": ["C5,4", "C4,5"]
}

def get_shortest_path(start, end, visited):
    if start == end:
        return [start]
    
    queue = deque([(start, [start])])
    seen = {start}
    
    while queue:
        current, path = queue.popleft()
        for next_pos in adjacency[current]:
            if next_pos not in seen and next_pos not in obstacles and next_pos not in visited:
                new_path = path + [next_pos]
                if next_pos == end:
                    return new_path
                seen.add(next_pos)
                queue.append((next_pos, new_path))
    return None

def find_path():
    # Order goals by proximity to current position
    current = initial_pos
    path = [current]
    visited_goals = set()
    if current in goals:
        visited_goals.add(current)
    
    # Define a specific order of goals that minimizes backtracking
    goal_order = ['C3,3', 'C3,2', 'C4,2', 'C5,1', 'C5,3', 'C5,4']
    
    for goal in goal_order:
        if goal not in visited_goals:
            segment = get_shortest_path(current, goal, set())
            if segment:
                path.extend(segment[1:])  # Exclude the start position to avoid duplication
                current = goal
                visited_goals.add(goal)
    
    return path

# Find and print the path
path = find_path()
print(f"<<<{json.dumps(path)}>>>")