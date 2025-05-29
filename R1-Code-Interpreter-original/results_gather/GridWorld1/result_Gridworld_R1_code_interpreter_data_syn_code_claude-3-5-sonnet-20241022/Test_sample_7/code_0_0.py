from collections import deque
import json

def find_path(start, end, adjacency, obstacles):
    if start == end:
        return [start]
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        for next_pos in adjacency[current]:
            if next_pos not in visited and next_pos not in obstacles:
                if next_pos == end:
                    return path + [next_pos]
                visited.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def find_complete_path():
    # Initialize data
    start = "C4,2"
    goals = ["C2,4", "C1,4", "C2,2", "C4,3"]
    obstacles = ["C3,3", "C3,1", "C4,1", "C2,3", "C1,3"]
    
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
    
    # Try different goal orders to find a good path
    current_pos = start
    remaining_goals = goals.copy()
    final_path = [start]
    
    # First, try to reach C2,2 as it's more central
    if "C2,2" in remaining_goals:
        next_goal = "C2,2"
        remaining_goals.remove(next_goal)
        path = find_path(current_pos, next_goal, adjacency, obstacles)
        if path:
            final_path.extend(path[1:])
            current_pos = next_goal
    
    # Then try to reach C4,3 as it's closer to the start
    if "C4,3" in remaining_goals:
        next_goal = "C4,3"
        remaining_goals.remove(next_goal)
        path = find_path(current_pos, next_goal, adjacency, obstacles)
        if path:
            final_path.extend(path[1:])
            current_pos = next_goal
    
    # Handle remaining goals
    while remaining_goals:
        shortest_path = None
        best_goal = None
        
        for goal in remaining_goals:
            path = find_path(current_pos, goal, adjacency, obstacles)
            if path and (shortest_path is None or len(path) < len(shortest_path)):
                shortest_path = path
                best_goal = goal
        
        if shortest_path:
            final_path.extend(shortest_path[1:])
            current_pos = best_goal
            remaining_goals.remove(best_goal)
        else:
            break
    
    print(json.dumps(final_path))

find_complete_path()