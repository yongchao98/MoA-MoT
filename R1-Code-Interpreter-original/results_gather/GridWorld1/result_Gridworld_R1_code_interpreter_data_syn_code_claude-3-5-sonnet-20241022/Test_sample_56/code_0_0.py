import json
from collections import deque

def find_path_visiting_all_goals(adjacency, start, goals, obstacles):
    def get_path_to_goals(current_pos, visited_goals):
        queue = deque([(current_pos, [current_pos], set(visited_goals))])
        seen = set()
        
        while queue:
            pos, path, current_visited_goals = queue.popleft()
            state = (pos, frozenset(current_visited_goals))
            
            if state in seen:
                continue
            seen.add(state)
            
            # Update visited goals
            if pos in goals:
                current_visited_goals = set(current_visited_goals)
                current_visited_goals.add(pos)
                
                # If all goals are visited, return the path
                if len(current_visited_goals) == len(goals):
                    return path
            
            # Try all possible moves
            for next_pos in adjacency[pos]:
                if next_pos not in obstacles:
                    queue.append((next_pos, path + [next_pos], current_visited_goals))
        
        return None

    # Find the path
    path = get_path_to_goals(start, set())
    return path

# Define the problem
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

start = "C3,4"
goals = ["C3,2", "C1,3", "C1,1", "C2,2"]
obstacles = ["C4,3", "C1,4", "C4,1", "C2,4", "C1,2"]

# Find the path
path = find_path_visiting_all_goals(adjacency, start, goals, obstacles)

# Format the output
output = json.dumps(path)
print(f"<<<{output}>>>")