import json
from collections import deque

def find_path(start, target, adjacency, obstacles):
    queue = deque([(start, [start])])
    seen = {start}
    
    while queue:
        current, path = queue.popleft()
        if current == target:
            return path
            
        for next_pos in adjacency[current]:
            if next_pos not in seen and next_pos not in obstacles:
                seen.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def solve_grid(start, goals, adjacency, obstacles):
    # Define the sequence of goals in a specific order
    goal_sequence = ["C5,4", "C5,3", "C3,3", "C3,2", "C4,2", "C5,1"]
    
    current = start
    final_path = [start]
    visited = {start}
    
    for goal in goal_sequence:
        # For each goal, find a path that avoids obstacles and previously visited cells
        # (except the current position)
        temp_obstacles = obstacles | (visited - {current})
        path = find_path(current, goal, adjacency, temp_obstacles)
        
        if not path:
            # If no path found, try without avoiding visited cells
            path = find_path(current, goal, adjacency, obstacles)
            
        if not path:
            return None
            
        # Add the new path (excluding the start position)
        final_path.extend(path[1:])
        visited.update(path[1:])
        current = goal
    
    return final_path

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

solution = solve_grid(start, goals, adjacency, obstacles)
if solution:
    print(f"<<<{json.dumps(solution)}>>>")
else:
    print("No valid path found")