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

def find_path_to_goals():
    def bfs(start, unvisited_goals):
        queue = deque([(start, [start], frozenset([start] if start in goals else []))])
        seen = {(start, frozenset([start] if start in goals else []))}
        
        while queue:
            pos, path, visited_goals = queue.popleft()
            
            # If we've visited all goals, return this path
            if all(g in visited_goals for g in goals):
                return path
                
            for next_pos in adjacency[pos]:
                if next_pos in obstacles:
                    continue
                    
                new_goals = set(visited_goals)
                if next_pos in goals:
                    new_goals.add(next_pos)
                new_goals_frozen = frozenset(new_goals)
                
                state = (next_pos, new_goals_frozen)
                if state not in seen:
                    seen.add(state)
                    queue.append((next_pos, path + [next_pos], new_goals_frozen))
        return None

    # Start the search from the initial position
    path = bfs(initial_pos, set(goals))
    
    # Remove any consecutive duplicates
    if path:
        optimized = [path[0]]
        for pos in path[1:]:
            if pos != optimized[-1]:
                optimized.append(pos)
        return optimized
    return None

# Find and print the path
path = find_path_to_goals()
print(f"<<<{json.dumps(path)}>>>")