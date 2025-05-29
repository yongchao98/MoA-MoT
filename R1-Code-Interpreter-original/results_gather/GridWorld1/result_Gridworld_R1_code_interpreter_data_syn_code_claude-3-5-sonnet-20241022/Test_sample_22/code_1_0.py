from collections import deque
import json

def find_path(start, end, adjacency, obstacles):
    queue = deque([(start, [start])])
    visited = set([start])
    
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
    current = start
    remaining_goals = set(goals)
    final_path = [start]
    
    while remaining_goals:
        shortest_path = None
        nearest_goal = None
        shortest_length = float('inf')
        
        for goal in remaining_goals:
            path = find_path(current, goal, adjacency, obstacles)
            if path and len(path) < shortest_length:
                shortest_length = len(path)
                shortest_path = path
                nearest_goal = goal
        
        if not shortest_path:
            return None
            
        final_path.extend(shortest_path[1:])
        current = nearest_goal
        remaining_goals.remove(nearest_goal)
    
    return final_path

# Problem setup
start = "C3,3"
goals = ["C4,1", "C4,2", "C2,2", "C1,2"]
obstacles = ["C1,5", "C3,4", "C1,1", "C4,3", "C4,5"]

# Complete adjacency dictionary
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

# Find the solution
solution = find_complete_path(start, goals, adjacency, obstacles)

# Output the solution
print(f"<<<{json.dumps(solution)}>>>")