import json
from collections import deque

def find_direct_path(start, goals, obstacles, adjacency):
    def bfs_with_constraints(start, target, must_visit=None, avoid=None):
        if must_visit is None:
            must_visit = set()
        if avoid is None:
            avoid = set()
            
        queue = deque([(start, [start], set())])
        visited = {start}
        
        while queue:
            current, path, visited_goals = queue.popleft()
            
            # Update visited goals
            if current in goals:
                visited_goals = visited_goals | {current}
                
            # Check if we've reached the target and visited all required points
            if current == target and must_visit.issubset(visited_goals):
                return path
                
            for next_pos in adjacency[current]:
                if (next_pos not in visited and 
                    next_pos not in obstacles and 
                    next_pos not in avoid):
                    visited.add(next_pos)
                    queue.append((next_pos, path + [next_pos], visited_goals))
        return None

    # Try different goal orderings
    possible_sequences = [
        ["C3,5", "C4,5", "C4,4", "C3,2", "C3,1", "C2,4", "C1,4"],  # Up then across
        ["C2,4", "C1,4", "C3,5", "C4,5", "C4,4", "C3,2", "C3,1"],  # First to C1,4
    ]
    
    best_path = None
    min_length = float('inf')
    
    for sequence in possible_sequences:
        current = start
        current_path = [start]
        path_valid = True
        visited = {start}
        
        for target in sequence:
            if target not in visited:  # Only move to unvisited positions
                path = bfs_with_constraints(current, target, avoid=visited-set(goals))
                if path is None:
                    path_valid = False
                    break
                current_path.extend(path[1:])  # Add new positions
                visited.update(path)
                current = target
        
        if path_valid and all(goal in current_path for goal in goals):
            if len(current_path) < min_length:
                min_length = len(current_path)
                best_path = current_path

    return best_path

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
path = find_direct_path(start, goals, obstacles, adjacency)

# Output the result
print(f"<<<{json.dumps(path)}>>>")