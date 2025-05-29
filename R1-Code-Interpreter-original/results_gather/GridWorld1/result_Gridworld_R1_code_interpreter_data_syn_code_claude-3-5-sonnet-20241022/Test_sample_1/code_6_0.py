import json
from collections import deque

def find_valid_path(start, goals, obstacles, adjacency):
    def bfs(start, target, visited_set=None):
        if visited_set is None:
            visited_set = set()
        queue = deque([(start, [start])])
        visited = {start}
        
        while queue:
            current, path = queue.popleft()
            if current == target:
                return path
            
            for next_pos in adjacency[current]:
                if (next_pos not in visited and 
                    next_pos not in obstacles):
                    visited.add(next_pos)
                    queue.append((next_pos, path + [next_pos]))
        return None

    # Start with closest goal first
    current = start
    final_path = [start]
    remaining_goals = set(goals)
    
    # First try to reach C1,4 since it's closest to start
    if "C1,4" in remaining_goals:
        path = bfs(current, "C1,4")
        if path:
            final_path.extend(path[1:])
            current = "C1,4"
            remaining_goals.remove("C1,4")
    
    # Then try to reach C4,5 and C4,4 since they're close to each other
    if "C4,5" in remaining_goals:
        path = bfs(current, "C4,5")
        if path:
            final_path.extend(path[1:])
            current = "C4,5"
            remaining_goals.remove("C4,5")
    
    if "C4,4" in remaining_goals:
        path = bfs(current, "C4,4")
        if path:
            final_path.extend(path[1:])
            current = "C4,4"
            remaining_goals.remove("C4,4")
    
    # Finally reach C3,1
    if "C3,1" in remaining_goals:
        path = bfs(current, "C3,1")
        if path:
            final_path.extend(path[1:])
            current = "C3,1"
            remaining_goals.remove("C3,1")
    
    # Check if all goals were reached
    if not remaining_goals:
        return final_path
    return None

# Define the problem parameters
start = "C2,5"
goals = ["C3,1", "C4,4", "C4,5", "C1,4"]
obstacles = ["C3,4", "C1,5", "C1,2", "C2,2", "C4,1"]
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

# Find the path
path = find_valid_path(start, goals, obstacles, adjacency)

# Output the result
print(f"<<<{json.dumps(path)}>>>")