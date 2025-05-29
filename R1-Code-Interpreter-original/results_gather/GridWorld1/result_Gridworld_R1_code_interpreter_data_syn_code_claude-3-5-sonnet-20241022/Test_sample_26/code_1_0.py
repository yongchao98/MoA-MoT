from collections import deque
import json

def find_shortest_path(adjacency, start, end, obstacles):
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

def find_complete_path(adjacency, start, goals, obstacles):
    current_pos = start
    remaining_goals = set(goals)
    path = [start]
    
    while remaining_goals:
        min_distance = float('inf')
        next_goal = None
        next_path = None
        
        for goal in remaining_goals:
            current_path = find_shortest_path(adjacency, current_pos, goal, obstacles)
            if current_path and len(current_path) < min_distance:
                min_distance = len(current_path)
                next_goal = goal
                next_path = current_path
        
        if next_goal is None:
            return None
            
        path.extend(next_path[1:])
        current_pos = next_goal
        remaining_goals.remove(next_goal)
    
    return path

# Problem setup
start = "C2,6"
goals = ['C6,3', 'C4,2', 'C2,5', 'C4,3', 'C4,5', 'C1,1', 'C3,3']
obstacles = ['C5,6', 'C4,6', 'C1,4', 'C1,3', 'C6,2', 'C3,4', 'C6,1', 'C3,2', 'C4,1', 'C6,6']

# Adjacency dictionary
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
    "C5,1": ["C5,2", "C4,1", "C6,1"],
    "C5,2": ["C5,1", "C5,3", "C4,2", "C6,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3", "C6,3"],
    "C5,4": ["C5,3", "C5,5", "C4,4", "C6,4"],
    "C5,5": ["C5,4", "C5,6", "C4,5", "C6,5"],
    "C5,6": ["C5,5", "C4,6", "C6,6"],
    "C6,1": ["C6,2", "C5,1"],
    "C6,2": ["C6,1", "C6,3", "C5,2"],
    "C6,3": ["C6,2", "C6,4", "C5,3"],
    "C6,4": ["C6,3", "C6,5", "C5,4"],
    "C6,5": ["C6,4", "C6,6", "C5,5"],
    "C6,6": ["C6,5", "C5,6"]
}

# Find the solution path
solution_path = find_complete_path(adjacency, start, goals, obstacles)

# Format the output
if solution_path:
    print(f"<<<{json.dumps(solution_path)}>>>")
else:
    print("No valid path found")