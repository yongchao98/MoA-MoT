import json
from collections import deque

def find_path(start, target, adjacency, obstacles, visited):
    queue = deque([(start, [start])])
    seen = {start}
    
    while queue:
        current, path = queue.popleft()
        if current == target:
            return path
            
        for next_pos in adjacency[current]:
            if (next_pos not in seen and 
                next_pos not in obstacles and 
                next_pos not in visited):
                seen.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def find_path_through_goals(start, goals, adjacency, obstacles):
    # Order goals based on position to minimize path crossings
    ordered_goals = []
    current = start
    remaining_goals = set(goals)
    path = [start]
    visited = set([start])
    
    # First visit C5,4 and C5,3
    if "C5,4" in remaining_goals:
        ordered_goals.append("C5,4")
        remaining_goals.remove("C5,4")
    if "C5,3" in remaining_goals:
        ordered_goals.append("C5,3")
        remaining_goals.remove("C5,3")
    
    # Then visit C3,3 and C3,2
    if "C3,3" in remaining_goals:
        ordered_goals.append("C3,3")
        remaining_goals.remove("C3,3")
    if "C3,2" in remaining_goals:
        ordered_goals.append("C3,2")
        remaining_goals.remove("C3,2")
    
    # Finally visit C4,2 and C5,1
    if "C4,2" in remaining_goals:
        ordered_goals.append("C4,2")
        remaining_goals.remove("C4,2")
    if "C5,1" in remaining_goals:
        ordered_goals.append("C5,1")
        remaining_goals.remove("C5,1")
    
    # Follow the predetermined order
    for goal in ordered_goals:
        temp_path = find_path(current, goal, adjacency, obstacles, visited - {current})
        if not temp_path:
            return None
        
        # Add new positions to path and visited set
        path.extend(temp_path[1:])
        visited.update(temp_path[1:])
        current = goal
    
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

solution = find_path_through_goals(start, goals, adjacency, obstacles)
if solution:
    print(f"<<<{json.dumps(solution)}>>>")
else:
    print("No valid path found")