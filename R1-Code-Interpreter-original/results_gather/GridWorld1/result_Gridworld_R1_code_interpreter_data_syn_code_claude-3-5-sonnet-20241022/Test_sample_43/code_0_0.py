import json
from collections import deque

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
    path = [start]
    
    while remaining_goals:
        # Find nearest goal
        min_path = None
        chosen_goal = None
        
        for goal in remaining_goals:
            if current == goal:
                chosen_goal = goal
                min_path = [current]
                break
                
            temp_path = find_path(current, goal, adjacency, obstacles)
            if temp_path and (min_path is None or len(temp_path) < len(min_path)):
                min_path = temp_path
                chosen_goal = goal
        
        if min_path is None:
            return None
            
        # Add path to goal (excluding start position if not first goal)
        path.extend(min_path[1:])
        remaining_goals.remove(chosen_goal)
        current = chosen_goal
    
    return path

# Define the problem
start = "C4,4"
goals = ["C2,4", "C1,5", "C3,3", "C4,5"]
obstacles = ["C3,1", "C4,1", "C2,2", "C1,2", "C3,5"]

# Define adjacency (your provided adjacency dict goes here)
adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    # ... (rest of the adjacency dict)
}

# Find the solution
solution = find_complete_path(start, goals, adjacency, obstacles)

# Format and print the solution
if solution:
    print(f"<<<{json.dumps(solution)}>>>")
else:
    print("No solution found")