import json
from collections import deque

def find_sequential_path(start, goals, obstacles, adjacency):
    def bfs(start, target):
        queue = deque([(start, [start])])
        visited = {start}
        
        while queue:
            current, path = queue.popleft()
            if current == target:
                return path
            
            for next_pos in adjacency[current]:
                if next_pos not in visited and next_pos not in obstacles:
                    visited.add(next_pos)
                    queue.append((next_pos, path + [next_pos]))
        return None

    # Define a specific sequence of goals that makes geometric sense
    goal_sequence = [
        "C2,4",  # First move towards C1,4
        "C1,4",  # First goal
        "C3,5",  # Move towards C4,5
        "C4,5",  # Second goal
        "C4,4",  # Third goal
        "C3,2",  # Move towards C3,1
        "C3,1"   # Fourth goal
    ]

    current = start
    final_path = [start]
    
    # Follow the sequence
    for target in goal_sequence:
        if current != target:  # Only find path if we're not already at the target
            path = bfs(current, target)
            if path is None:
                return None  # If any segment fails, the whole path is invalid
            final_path.extend(path[1:])  # Add new positions (excluding current position)
            current = target

    # Verify that all goals are in the path
    if not all(goal in final_path for goal in goals):
        return None

    return final_path

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
path = find_sequential_path(start, goals, obstacles, adjacency)

# Output the result
print(f"<<<{json.dumps(path)}>>>")