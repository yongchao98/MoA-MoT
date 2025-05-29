import json
from collections import deque
import copy

# Initialize the problem
initial_pos = "C5,3"
goals = ['C5,4', 'C4,5', 'C1,3', 'C3,4', 'C3,3', 'C1,2']
obstacles = ['C1,4', 'C3,5', 'C4,3', 'C5,2', 'C5,1', 'C1,5', 'C2,1', 'C5,5']

# Adjacency list
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

def find_path():
    # Queue will store: (current_position, path_so_far, visited_goals)
    queue = deque([(initial_pos, [initial_pos], set([initial_pos]) if initial_pos in goals else set())])
    seen_states = set()  # To avoid cycles: (position, frozenset(visited_goals))
    
    while queue:
        current, path, visited_goals = queue.popleft()
        
        # If we've visited all goals, return the path
        if all(goal in visited_goals for goal in goals):
            return path
            
        # Get valid next positions
        for next_pos in adjacency[current]:
            if next_pos in obstacles:
                continue
                
            # Calculate new visited goals
            new_visited_goals = visited_goals | ({next_pos} if next_pos in goals else set())
            
            # Create state key to check if we've seen this state before
            state_key = (next_pos, frozenset(new_visited_goals))
            if state_key in seen_states:
                continue
                
            seen_states.add(state_key)
            new_path = path + [next_pos]
            queue.append((next_pos, new_path, new_visited_goals))
    
    return None

# Find and print solution
solution = find_path()
if solution:
    print(f"<<<{json.dumps(solution)}>>>")
else:
    print("No solution found")