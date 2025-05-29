import json

initial_state = {
    "box1": "C1,4",
    "box2": "C2,8",
    "box3": "C2,2",
    "box4": "C2,1"
}

goal_state = {
    "box1": "C2,8",
    "box2": "C2,2",
    "box3": "C1,2",
    "box4": "C2,5"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"],
    "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"],
    "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C1,6", "C2,5"],
    "C1,6": ["C1,5", "C1,7", "C2,6"],
    "C1,7": ["C1,6", "C1,8", "C2,7"],
    "C1,8": ["C1,7", "C2,8"],
    "C2,1": ["C2,2", "C1,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2"],
    "C2,3": ["C2,2", "C2,4", "C1,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4"],
    "C2,5": ["C2,4", "C2,6", "C1,5"],
    "C2,6": ["C2,5", "C2,7", "C1,6"],
    "C2,7": ["C2,6", "C2,8", "C1,7"],
    "C2,8": ["C2,7", "C1,8"]
}

def move_box_towards_goal(current_state, box):
    states = []
    current = current_state.copy()
    
    while current[box] != goal_state[box]:
        current_pos = current[box]
        target_pos = goal_state[box]
        
        # Get row and column of current and target positions
        curr_row = int(current_pos.split(',')[0][1])
        curr_col = int(current_pos.split(',')[1])
        target_row = int(target_pos.split(',')[0][1])
        target_col = int(target_pos.split(',')[1])
        
        # Try to move towards target
        occupied = set(current.values())
        next_pos = None
        
        # Try to move horizontally first
        if curr_col < target_col:
            next_candidate = f"C{curr_row},{curr_col+1}"
            if next_candidate in adjacency[current_pos] and next_candidate not in occupied:
                next_pos = next_candidate
        elif curr_col > target_col:
            next_candidate = f"C{curr_row},{curr_col-1}"
            if next_candidate in adjacency[current_pos] and next_candidate not in occupied:
                next_pos = next_candidate
                
        # If can't move horizontally, try vertically
        if not next_pos:
            if curr_row < target_row:
                next_candidate = f"C{curr_row+1},{curr_col}"
                if next_candidate in adjacency[current_pos] and next_candidate not in occupied:
                    next_pos = next_candidate
            elif curr_row > target_row:
                next_candidate = f"C{curr_row-1},{curr_col}"
                if next_candidate in adjacency[current_pos] and next_candidate not in occupied:
                    next_pos = next_candidate
        
        if next_pos:
            current[box] = next_pos
            states.append(current.copy())
        else:
            break
            
    return states

def solve():
    solution = [initial_state]
    current_state = initial_state.copy()
    
    # Move boxes in sequence: box1 -> box2 -> box3 -> box4
    for box in ["box1", "box2", "box3", "box4"]:
        if current_state[box] != goal_state[box]:
            next_states = move_box_towards_goal(current_state, box)
            if next_states:
                solution.extend(next_states)
                current_state = next_states[-1]
    
    return solution

# Generate and print solution
solution = solve()
print(json.dumps(solution))