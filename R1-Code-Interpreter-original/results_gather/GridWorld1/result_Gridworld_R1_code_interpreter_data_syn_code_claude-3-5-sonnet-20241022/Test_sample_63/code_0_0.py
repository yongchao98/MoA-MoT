from collections import deque
import json

def find_path_to_all_goals():
    start = "C2,3"
    goals = {"C4,4", "C3,2", "C4,3", "C4,1"}
    obstacles = {"C2,1", "C3,4", "C1,3", "C1,2", "C3,1"}
    
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
    
    # State: (current_position, frozenset of visited_goals, path_so_far)
    queue = deque([(start, frozenset(), [start])])
    seen = set()
    
    while queue:
        pos, visited_goals, path = queue.popleft()
        state = (pos, visited_goals)
        
        if state in seen:
            continue
        seen.add(state)
        
        # Update visited goals
        if pos in goals:
            visited_goals = frozenset(list(visited_goals) + [pos])
            
        # Check if all goals are visited
        if len(visited_goals) == len(goals):
            return path
            
        # Try all possible moves
        for next_pos in adjacency[pos]:
            if next_pos not in obstacles:
                queue.append((next_pos, visited_goals, path + [next_pos]))
    
    return None

# Find and print the path
path = find_path_to_all_goals()
print(f"<<<{json.dumps(path)}>>>")