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

def find_optimal_path():
    def bfs():
        start_state = (initial_pos, frozenset())
        queue = deque([(start_state, [initial_pos])])
        seen = {start_state}
        
        while queue:
            (current_pos, visited_goals), path = queue.popleft()
            
            # Update visited goals
            new_goals = set(visited_goals)
            if current_pos in goals:
                new_goals.add(current_pos)
            
            # Check if we've found a solution
            if len(new_goals) == len(goals):
                return path
            
            # Try each possible move
            for next_pos in adjacency[current_pos]:
                if next_pos in obstacles:
                    continue
                
                new_state = (next_pos, frozenset(new_goals))
                if new_state not in seen:
                    seen.add(new_state)
                    queue.append((new_state, path + [next_pos]))
        
        return None

    # Find the shortest path visiting all goals
    path = bfs()
    
    if not path:
        return None
    
    # Remove any unnecessary moves
    final_path = [path[0]]
    for i in range(1, len(path)):
        if path[i] != final_path[-1]:
            final_path.append(path[i])
    
    return final_path

# Find and print the path
path = find_optimal_path()
print(f"<<<{json.dumps(path)}>>>")