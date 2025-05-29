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

def is_valid_move(state, box, new_pos):
    # Check if the new position is occupied by any other box
    for other_box, pos in state.items():
        if other_box != box and pos == new_pos:
            return False
    return True

def solve():
    solution = []
    current = initial_state.copy()
    solution.append(current.copy())
    
    # Move sequence:
    # 1. First move box3 to C1,2 (to clear space for other moves)
    current = current.copy()
    current["box3"] = "C1,2"
    solution.append(current.copy())
    
    # 2. Move box1 to C2,8
    moves_box1 = ["C1,5", "C1,6", "C1,7", "C1,8", "C2,8"]
    for pos in moves_box1:
        if is_valid_move(current, "box1", pos):
            current = current.copy()
            current["box1"] = pos
            solution.append(current.copy())
    
    # 3. Move box2 to C2,2
    moves_box2 = ["C2,7", "C2,6", "C2,5", "C2,4", "C2,3", "C2,2"]
    for pos in moves_box2:
        if is_valid_move(current, "box2", pos):
            current = current.copy()
            current["box2"] = pos
            solution.append(current.copy())
    
    # 4. Move box4 to C2,5
    moves_box4 = ["C2,2", "C2,3", "C2,4", "C2,5"]
    for pos in moves_box4:
        if is_valid_move(current, "box4", pos):
            current = current.copy()
            current["box4"] = pos
            solution.append(current.copy())
    
    # Verify final state matches goal state
    if current == goal_state:
        return solution
    return None

# Generate and print solution
solution = solve()
if solution:
    print(json.dumps(solution))
else:
    print("No valid solution found")