def find_valid_path(start, goals, adjacency, obstacles):
    goals_set = set(goals)
    best_path = None
    
    def dfs(current, path, visited_goals):
        nonlocal best_path
        
        # If we already found a path and current path is longer, skip
        if best_path is not None and len(path) >= len(best_path):
            return
            
        # Update visited goals
        current_visited = visited_goals | ({current} if current in goals_set else set())
        
        # If we've visited all goals, we found a valid path
        if len(current_visited) == len(goals_set):
            best_path = path
            return
            
        # Try each possible move
        for next_pos in adjacency[current]:
            if next_pos not in obstacles and next_pos not in path:
                dfs(next_pos, path + [next_pos], current_visited)
    
    # Start DFS from initial position
    dfs(start, [start], set())
    
    return best_path

# Problem setup
start = "C4,3"
goals = ["C3,4", "C3,3", "C2,4", "C2,1"]
obstacles = ["C2,2", "C4,4", "C2,3", "C1,2", "C4,2"]
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

result = find_valid_path(start, goals, adjacency, obstacles)
if result:
    print(result)