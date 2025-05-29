import json
from collections import deque

# Initial and goal states
initial_state = {
    "box1": "C1,7", "box2": "C2,8", "box3": "C4,5",
    "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"
}
goal_state = {
    "box1": "C3,6", "box2": "C4,1", "box3": "C2,7",
    "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"
}

# Adjacency dictionary
adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"], "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"], "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1", "C3,1"], "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"], "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"], "C2,6": ["C2,5", "C2,7", "C1,6", "C3,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7", "C3,7"], "C2,8": ["C2,7", "C1,8", "C3,8"],
    "C3,1": ["C3,2", "C2,1", "C4,1"], "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"], "C3,6": ["C3,5", "C3,7", "C2,6", "C4,6"],
    "C3,7": ["C3,6", "C3,8", "C2,7", "C4,7"], "C3,8": ["C3,7", "C2,8", "C4,8"],
    "C4,1": ["C4,2", "C3,1"], "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"], "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C4,5": ["C4,4", "C4,6", "C3,5"], "C4,6": ["C4,5", "C4,7", "C3,6"],
    "C4,7": ["C4,6", "C4,8", "C3,7"], "C4,8": ["C4,7", "C3,8"]
}

def find_path(start, end, blocked):
    if start == end:
        return []
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        pos, path = queue.popleft()
        for next_pos in adjacency[pos]:
            if next_pos not in visited and next_pos not in blocked:
                if next_pos == end:
                    return path + [next_pos]
                visited.add(next_pos)
                queue.append((next_pos, path + [next_pos]))
    return None

def move_box(state, box, path):
    if not path or len(path) < 2:
        return None
    result = []
    current_state = state.copy()
    
    for i in range(1, len(path)):
        new_state = current_state.copy()
        new_state[box] = path[i]
        result.append(new_state)
        current_state = new_state
    
    return result

def solve():
    # Movement order based on final positions
    box_order = ["box1", "box6", "box3", "box4", "box5", "box2"]
    current_state = initial_state.copy()
    path = [current_state]
    
    for box in box_order:
        target = goal_state[box]
        while current_state[box] != target:
            blocked = set(current_state.values()) - {current_state[box]}
            box_path = find_path(current_state[box], target, blocked)
            
            if not box_path:
                # Try to move other boxes out of the way
                for other_box in current_state:
                    if other_box != box and other_box not in [b for b in box_order[:box_order.index(box)]]:
                        other_pos = current_state[other_box]
                        for alt_pos in adjacency[other_pos]:
                            if alt_pos not in blocked and alt_pos != target:
                                new_state = current_state.copy()
                                new_state[other_box] = alt_pos
                                path.append(new_state)
                                current_state = new_state
                                blocked = set(current_state.values()) - {current_state[box]}
                                box_path = find_path(current_state[box], target, blocked)
                                if box_path:
                                    break
                        if box_path:
                            break
                
                if not box_path:
                    return None
            
            moves = move_box(current_state, box, box_path)
            if not moves:
                return None
            
            path.extend(moves)
            current_state = moves[-1]
            
            if len(path) > 100:  # Prevent infinite loops
                return None
    
    return path

# Find and print solution
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")