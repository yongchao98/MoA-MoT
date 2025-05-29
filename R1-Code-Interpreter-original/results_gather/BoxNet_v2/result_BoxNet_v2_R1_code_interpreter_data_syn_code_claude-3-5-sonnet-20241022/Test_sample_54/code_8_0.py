import json

def is_valid_state(state):
    # Check if any two boxes share the same position
    positions = list(state.values())
    return len(positions) == len(set(positions))

def direct_solve():
    path = []
    current_state = {
        "box1": "C1,4",
        "box2": "C2,8",
        "box3": "C2,2",
        "box4": "C2,1"
    }
    path.append(current_state.copy())

    # Move box1 to C2,8 (first move it to C1,8, then to C2,8)
    moves_box1 = ["C1,5", "C1,6", "C1,7", "C1,8"]
    for move in moves_box1:
        new_state = current_state.copy()
        new_state["box1"] = move
        path.append(new_state.copy())
        current_state = new_state

    # Move box1 from C1,8 to C2,8 (only after box2 moves)
    current_state = path[-1].copy()
    new_state = current_state.copy()
    new_state["box2"] = "C2,7"  # Move box2 first
    path.append(new_state.copy())
    current_state = new_state

    new_state = current_state.copy()
    new_state["box1"] = "C2,8"  # Then move box1
    path.append(new_state.copy())
    current_state = new_state

    # Move box2 to C2,2
    moves_box2 = ["C2,6", "C2,5", "C2,4", "C2,3", "C2,2"]
    for move in moves_box2:
        new_state = current_state.copy()
        new_state["box2"] = move
        path.append(new_state.copy())
        current_state = new_state

    # Move box3 to C1,2
    new_state = current_state.copy()
    new_state["box3"] = "C1,2"
    path.append(new_state.copy())
    current_state = new_state

    # Move box4 to C2,5
    moves_box4 = ["C2,2", "C2,3", "C2,4", "C2,5"]
    for move in moves_box4:
        new_state = current_state.copy()
        new_state["box4"] = move
        path.append(new_state.copy())
        current_state = new_state

    # Verify all states are valid
    valid_path = []
    for state in path:
        if is_valid_state(state):
            valid_path.append(state)

    return valid_path

# Get solution and print
solution = direct_solve()
print(json.dumps(solution))