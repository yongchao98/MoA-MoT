from collections import deque
import json

# Initialize the problem
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C2,4"],
    "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C1,4", "C3,4"],
    "C3,1": ["C3,2", "C2,1", "C4,1"],
    "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C3,4": ["C3,3", "C2,4", "C4,4"],
    "C4,1": ["C4,2", "C3,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"],
    "C4,4": ["C4,3", "C3,4"]
}

obstacles = ['C2,3', 'C3,3', 'C3,2', 'C2,4', 'C3,4']
goals = set(['C1,2', 'C4,1', 'C1,3', 'C4,4'])
start = 'C1,1'

def find_path():
    # Queue will store (current_pos, path_so_far, goals_visited)
    queue = deque([(start, [start], set())])
    visited = set()  # Store (pos, goals_visited) pairs
    
    while queue:
        current, path, goals_reached = queue.popleft()
        state = (current, tuple(sorted(goals_reached)))
        
        if state in visited:
            continue
        visited.add(state)
        
        # Update goals reached
        if current in goals:
            goals_reached = goals_reached | {current}
            
        # Check if all goals are reached
        if len(goals_reached) == len(goals):
            return path
            
        # Try all possible moves
        for next_pos in adjacency[current]:
            if next_pos not in obstacles:
                new_path = path + [next_pos]
                queue.append((next_pos, new_path, goals_reached))
    
    return None

# Find and print the solution
solution = find_path()
print(f"<<<{json.dumps(solution)}>>>")