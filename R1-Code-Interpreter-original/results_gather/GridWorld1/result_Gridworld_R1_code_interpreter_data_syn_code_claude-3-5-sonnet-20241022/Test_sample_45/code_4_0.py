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

def bfs_path(start, target, visited_goals):
    queue = deque([(start, [start], visited_goals)])
    seen = {(start, frozenset(visited_goals))}
    
    while queue:
        current, path, current_goals = queue.popleft()
        
        for next_pos in adjacency[current]:
            if next_pos in obstacles:
                continue
                
            new_goals = set(current_goals)
            if next_pos in goals:
                new_goals.add(next_pos)
            
            state = (next_pos, frozenset(new_goals))
            if state not in seen:
                new_path = path + [next_pos]
                if len(new_goals) == len(goals):
                    return new_path
                    
                seen.add(state)
                queue.append((next_pos, new_path, new_goals))
    return None

def find_optimal_path():
    initial_goals = {initial_pos} if initial_pos in goals else set()
    path = bfs_path(initial_pos, None, initial_goals)
    
    # Remove consecutive duplicates if any exist
    if path:
        optimized_path = [path[0]]
        for pos in path[1:]:
            if pos != optimized_path[-1]:
                optimized_path.append(pos)
        return optimized_path
    return None

# Find and print the path
path = find_optimal_path()
print(f"<<<{json.dumps(path)}>>>")