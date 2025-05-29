import json
from collections import deque

def find_path_to_all_goals(start, goals, adjacency, obstacles):
    def is_valid_move(pos):
        return pos not in obstacles
    
    def get_path_to_goal(start, current_goals, visited):
        queue = deque([(start, [start], set())])
        seen = set()
        
        while queue:
            current, path, visited_goals = queue.popleft()
            state = (current, tuple(sorted(visited_goals)))
            
            if state in seen:
                continue
            seen.add(state)
            
            # Update visited goals
            if current in current_goals:
                visited_goals = visited_goals | {current}
            
            # Check if all goals are visited
            if len(visited_goals) == len(current_goals):
                return path
            
            # Try all possible moves
            for next_pos in adjacency[current]:
                if is_valid_move(next_pos) and next_pos not in visited:
                    new_path = path + [next_pos]
                    queue.append((next_pos, new_path, visited_goals))
        
        return None

    # Initialize
    start_pos = start
    goals_set = set(goals)
    path = find_path_to_goal(start_pos, goals_set, set())
    
    if path:
        return json.dumps(path)
    return None

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

result = find_path_to_all_goals(start, goals, adjacency, obstacles)
print(result)