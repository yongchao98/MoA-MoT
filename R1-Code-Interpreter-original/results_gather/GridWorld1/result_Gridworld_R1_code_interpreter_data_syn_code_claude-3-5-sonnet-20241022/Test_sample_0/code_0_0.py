import json
from collections import deque

def find_path_to_all_goals():
    # Initialize data
    start = "C2,3"
    goals = {"C2,1", "C2,4", "C3,3", "C4,2"}
    obstacles = {"C1,3", "C2,2", "C1,2", "C4,4", "C1,1"}
    
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

    def get_valid_neighbors(pos):
        return [n for n in adjacency[pos] if n not in obstacles]

    def find_path(start_pos, target_pos, visited):
        if start_pos == target_pos:
            return [start_pos]
        
        queue = deque([(start_pos, [start_pos])])
        seen = {start_pos}
        
        while queue:
            current, path = queue.popleft()
            for next_pos in get_valid_neighbors(current):
                if next_pos not in seen and next_pos not in visited:
                    if next_pos == target_pos:
                        return path + [next_pos]
                    seen.add(next_pos)
                    queue.append((next_pos, path + [next_pos]))
        return None

    # Find path visiting all goals
    current_pos = start
    remaining_goals = goals.copy()
    final_path = [start]
    visited = set(obstacles)

    while remaining_goals:
        # Find nearest goal
        min_path = None
        min_goal = None
        
        for goal in remaining_goals:
            path = find_path(current_pos, goal, visited)
            if path and (min_path is None or len(path) < len(min_path)):
                min_path = path
                min_goal = goal

        if min_path is None:
            return None  # No valid path found
            
        # Add path to final path (excluding start position if not first path)
        final_path.extend(min_path[1:])
        current_pos = min_goal
        remaining_goals.remove(min_goal)

    print(json.dumps(final_path))

find_path_to_all_goals()