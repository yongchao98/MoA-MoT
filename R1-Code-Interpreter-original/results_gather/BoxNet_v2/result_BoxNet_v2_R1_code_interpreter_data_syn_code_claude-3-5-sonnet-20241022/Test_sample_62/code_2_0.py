import json

initial_state = {
    "box1": "C5,5", "box2": "C1,3", "box3": "C3,3",
    "box4": "C4,5", "box5": "C4,4"
}

goal_state = {
    "box1": "C2,5", "box2": "C3,4", "box3": "C1,5",
    "box4": "C3,1", "box5": "C1,3"
}

adjacency = {
    "C1,1": ["C1,2", "C2,1"], "C1,2": ["C1,1", "C1,3", "C2,2"],
    "C1,3": ["C1,2", "C1,4", "C2,3"], "C1,4": ["C1,3", "C1,5", "C2,4"],
    "C1,5": ["C1,4", "C2,5"], "C2,1": ["C2,2", "C1,1", "C3,1"],
    "C2,2": ["C2,1", "C2,3", "C1,2", "C3,2"], "C2,3": ["C2,2", "C2,4", "C1,3", "C3,3"],
    "C2,4": ["C2,3", "C2,5", "C1,4", "C3,4"], "C2,5": ["C2,4", "C1,5", "C3,5"],
    "C3,1": ["C3,2", "C2,1", "C4,1"], "C3,2": ["C3,1", "C3,3", "C2,2", "C4,2"],
    "C3,3": ["C3,2", "C3,4", "C2,3", "C4,3"], "C3,4": ["C3,3", "C3,5", "C2,4", "C4,4"],
    "C3,5": ["C3,4", "C2,5", "C4,5"], "C4,1": ["C4,2", "C3,1", "C5,1"],
    "C4,2": ["C4,1", "C4,3", "C3,2", "C5,2"], "C4,3": ["C4,2", "C4,4", "C3,3", "C5,3"],
    "C4,4": ["C4,3", "C4,5", "C3,4", "C5,4"], "C4,5": ["C4,4", "C3,5", "C5,5"],
    "C5,1": ["C5,2", "C4,1"], "C5,2": ["C5,1", "C5,3", "C4,2"],
    "C5,3": ["C5,2", "C5,4", "C4,3"], "C5,4": ["C5,3", "C5,5", "C4,4"],
    "C5,5": ["C5,4", "C4,5"]
}

def move_box(current_pos, target_pos, state, box_name):
    result = []
    current = current_pos
    
    # Extract row and column from positions
    curr_row, curr_col = map(int, current.replace('C', '').split(','))
    target_row, target_col = map(int, target_pos.replace('C', '').split(','))
    
    # Move row-wise first, then column-wise
    while curr_row != target_row or curr_col != target_col:
        new_state = dict(state)
        
        # Move row-wise
        if curr_row < target_row:
            curr_row += 1
        elif curr_row > target_row:
            curr_row -= 1
        # Move column-wise
        elif curr_col < target_col:
            curr_col += 1
        elif curr_col > target_col:
            curr_col -= 1
            
        new_pos = f"C{curr_row},{curr_col}"
        if new_pos in adjacency[current]:
            new_state[box_name] = new_pos
            result.append(new_state)
            current = new_pos
            state = new_state
            
    return result

def solve():
    solution = [initial_state]
    current_state = dict(initial_state)
    
    # Move box1 to goal
    moves = move_box(current_state["box1"], goal_state["box1"], current_state, "box1")
    for state in moves:
        solution.append(state)
        current_state = state
    
    # Move box2 to goal
    moves = move_box(current_state["box2"], goal_state["box2"], current_state, "box2")
    for state in moves:
        solution.append(state)
        current_state = state
        
    # Move box3 to goal
    moves = move_box(current_state["box3"], goal_state["box3"], current_state, "box3")
    for state in moves:
        solution.append(state)
        current_state = state
        
    # Move box4 to goal
    moves = move_box(current_state["box4"], goal_state["box4"], current_state, "box4")
    for state in moves:
        solution.append(state)
        current_state = state
        
    # Move box5 to goal
    moves = move_box(current_state["box5"], goal_state["box5"], current_state, "box5")
    for state in moves:
        solution.append(state)
        current_state = state
    
    return solution

# Generate and print solution
solution = solve()
print(json.dumps(solution))