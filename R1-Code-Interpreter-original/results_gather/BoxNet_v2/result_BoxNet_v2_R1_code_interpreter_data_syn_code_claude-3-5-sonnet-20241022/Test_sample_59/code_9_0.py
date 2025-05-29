import json

# Initial and goal states
initial_state = {
    "box1": "C1,7", "box2": "C2,8", "box3": "C4,5",
    "box4": "C2,4", "box5": "C4,7", "box6": "C3,6"
}
goal_state = {
    "box1": "C3,6", "box2": "C4,1", "box3": "C2,7",
    "box4": "C3,3", "box5": "C2,3", "box6": "C2,6"
}

# Simplified adjacency
adjacency = {
    "C1,7": ["C1,6", "C1,8", "C2,7"],
    "C2,8": ["C2,7", "C1,8", "C3,8"],
    "C4,5": ["C4,4", "C4,6", "C3,5"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"],
    "C4,7": ["C4,6", "C4,8", "C3,7"],
    "C3,6": ["C3,5", "C3,7", "C2,6", "C4,6"],
    "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,6": ["C2,5", "C2,7", "C1,6", "C3,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7", "C3,7"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"],
    "C4,1": ["C4,2", "C3,1"],
    "C3,8": ["C3,7", "C2,8", "C4,8"],
    "C4,6": ["C4,5", "C4,7", "C3,6"],
    "C3,7": ["C3,6", "C3,8", "C2,7", "C4,7"],
    "C4,8": ["C4,7", "C3,8"],
    "C4,2": ["C4,1", "C4,3", "C3,2"],
    "C4,3": ["C4,2", "C4,4", "C3,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4"],
    "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C3,6", "C2,5", "C4,5"],
    "C2,5": ["C2,4", "C2,6", "C1,5", "C3,5"]
}

def is_valid_move(state, box, new_pos):
    # Check if the new position is occupied
    if new_pos in state.values():
        return False
    return True

def try_move_box(state, box, path, visited):
    if len(path) > 20:  # Limit path length
        return None
        
    # If this box is at its goal, try next box
    if state[box] == goal_state[box]:
        return path
        
    current_pos = state[box]
    if current_pos not in adjacency:
        return None
        
    # Try each possible move for this box
    for next_pos in adjacency[current_pos]:
        if is_valid_move(state, box, next_pos):
            new_state = state.copy()
            new_state[box] = next_pos
            
            state_str = json.dumps(new_state, sort_keys=True)
            if state_str in visited:
                continue
                
            visited.add(state_str)
            new_path = path + [new_state]
            
            # If this move got us to the goal for this box, return the path
            if next_pos == goal_state[box]:
                return new_path
                
            # Otherwise, keep trying moves
            result = try_move_box(new_state, box, new_path, visited)
            if result:
                return result
                
            visited.remove(state_str)
    
    return None

def solve():
    current_state = initial_state.copy()
    full_path = [current_state]
    visited = {json.dumps(current_state, sort_keys=True)}
    
    # Move boxes one at a time to their goals
    for box in ["box1", "box6", "box4", "box5", "box3", "box2"]:
        path = try_move_box(current_state, box, full_path, visited)
        if not path:
            return None
        full_path = path
        current_state = path[-1]
    
    return full_path

# Find and print solution
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No solution found")