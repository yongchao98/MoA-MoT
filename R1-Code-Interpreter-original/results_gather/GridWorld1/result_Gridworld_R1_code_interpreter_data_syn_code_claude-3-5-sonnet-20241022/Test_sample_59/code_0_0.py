import json

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

def find_path(current, path, visited_goals, all_goals, obstacles, adjacency):
    # If we've visited all goals, we're done
    if all(goal in visited_goals for goal in all_goals):
        return path
    
    # Try each adjacent position
    for next_pos in adjacency[current]:
        # Skip if position is an obstacle
        if next_pos in obstacles:
            continue
            
        # Skip if we've already been here and we're not finding new goals
        if next_pos in path and next_pos not in all_goals:
            continue
            
        # Create new path and visited goals set
        new_path = path + [next_pos]
        new_visited_goals = visited_goals | ({next_pos} if next_pos in all_goals else set())
        
        # Recursively try this path
        result = find_path(next_pos, new_path, new_visited_goals, all_goals, obstacles, adjacency)
        if result:
            return result
            
    return None

# Find the solution
initial_visited_goals = {initial_pos} if initial_pos in goals else set()
solution = find_path(initial_pos, [initial_pos], initial_visited_goals, goals, obstacles, adjacency)

# Print the solution in the required format
if solution:
    print(f"<<<{json.dumps(solution)}>>>")
else:
    print("No solution found")