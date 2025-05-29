import json
from heapq import heappush, heappop

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

def manhattan_distance(pos1, pos2):
    x1, y1 = map(int, pos1.split(','))
    x2, y2 = map(int, pos2.split(','))
    return abs(x1 - x2) + abs(y1 - y2)

def find_path():
    # Priority queue for A* search
    # Format: (priority, path_length, current_pos, path, visited_goals)
    start_goals = frozenset([initial_pos]) if initial_pos in goals else frozenset()
    queue = [(0, 1, initial_pos, [initial_pos], start_goals)]
    seen = set()  # Format: (current_pos, visited_goals)
    
    while queue:
        _, path_len, current, path, visited_goals = heappop(queue)
        
        # Check if we've found a solution
        if len(visited_goals) == len(set(goals)):
            return path
            
        state = (current, visited_goals)
        if state in seen:
            continue
        seen.add(state)
        
        # Try all possible next moves
        for next_pos in adjacency[current]:
            if next_pos in obstacles:
                continue
                
            new_goals = visited_goals
            if next_pos in goals:
                new_goals = frozenset(visited_goals | {next_pos})
                
            if (next_pos, new_goals) in seen:
                continue
                
            # Calculate heuristic: minimum distance to any unvisited goal
            remaining_goals = set(goals) - set(new_goals)
            if remaining_goals:
                h = min(manhattan_distance(next_pos, g) for g in remaining_goals)
            else:
                h = 0
                
            new_path = path + [next_pos]
            priority = path_len + h
            
            heappush(queue, (priority, path_len + 1, next_pos, new_path, new_goals))
    
    return None

# Find and print the path
path = find_path()
print(f"<<<{json.dumps(path)}>>>")